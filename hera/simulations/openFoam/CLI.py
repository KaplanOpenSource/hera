import json
import logging
import glob
import os
import shutil

import pandas

from ...datalayer import datatypes
from ... import toolkitHome
from ...utils.jsonutils import loadJSON,compareJSONS
from ..hermesWorkflowToolkit import workflowsTypes

def stochasticLagrangian_create_dispersionFlow(arguments):
    logger = logging.getLogger("hera.bin")
    logger.execution(f"----- Start -----")
    logger.debug(f" makeDispersionFlow with arguments: {arguments}")

    if 'projectName' not in arguments:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"

        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']
    else:
        projectName = arguments.projectName

    logger.info(f"Adding dispersion flow to project {projectName}")
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    for flowid,flowdata in enumerate(configuration["flowFields"]["Flows"]):
        logger.info(f'Creating flow Dispersion {flowid}/{len(configuration["flowFields"]["Flows"])}')
        tk.stochasticLagrangian.createDispersionFlowField(flowData=flowdata)

    logger.execution(f"----- End -----")


def stochasticLagrangian_list_dispersionFlow(arguments):
    """
        Lists the dispersion flows and the differences between them
    Parameters
    ----------
    arguments

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin")
    logger.execution(f"----- Start -----")
    logger.debug(f" makeDispersionFlow with arguments: {arguments}")

    if 'projectName' not in arguments:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"

        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']
    else:
        projectName = arguments.projectName

    logger.info(f"Listing all dispersion flows in project {projectName}")
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)

    groupsList = pandas.DataFrame(dict(groups=[doc.desc['groupName'] for doc in tk.getSimulationDocuments(type=workflowsTypes.OF_FLOWDISPERSION.value)])).drop_duplicates()

    resList = []
    for group in groupsList:
        ttl = f"Group name {group}"
        print(ttl)
        print("-"*len(ttl))
        wfDict = dict([(doc.desc['name'],doc.desc['name']['flowParameters']) for doc in  tk.getSimulationDocuments(type=workflowsTypes.OF_FLOWDISPERSION.value,groupName=group)])
        res = compareJSONS(wfDict)
        if arguments.format == "pandas":
            print(res)
        else:
            resList.append(res.assign(groupName=group))

    if arguments.format != "pandas":
        allRes = pandas.concat(resList)

        if arguments.format == "latex":
            output = allRes.to_latex()
            ext = "tex"
        elif arguments.format == "csv":
            output = allRes.to_csv()
            ext = "csv"
        else:
            output = json.dumps(loadJSON(allRes.to_json()), indent=4)
            ext = "json"

        if len(allRes) == 0:
            print(f"Could not found any workflows to compare in project {projectName}")
        else:
            print(output)

            if arguments.file is not None:
                flName = arguments.file if "." in arguments.file else f"{arguments.file}.{ext}"

                with open(flName, "w") as outputFile:
                    outputFile.write(output)


def stochasticLagrangian_createDispersion(arguments):
    """
        Prepares the dispersion case:

        1. Create the dispersion case (and link to dispersionFlow)
        2. build and execute the dispersion workflow:
            - If the workflow on disk  mismatch the DB stop unless
                the --exportDB or the --updateDB flags exist.
            - If --exportDB flag exists, use the workflow in the DB.
        3. Add to the DB:
            - If the workflow mismatch the DB, update the DB if the --updateDB flag exists.
            - add the parameters of the flowcase to the DB.



    :param arguments:
             - dispersionCaseDirectory : the new case of the directory
             - dispersionFlow : the new case of the directory.
                                can be either a directory or a name of a dispersion flow field from the DB

    :return:
    """
    logger = logging.getLogger("hera.bin.stochasticLagrangian.createDispersionFlow")
    logger.execution(f"----- Start -----")
    logger.debug(f"  Got arguments: {arguments}")

    if (arguments.updateDB and arguments.exportFromDB):
        err = "Cannot use both --updateDB and --exportFromDB"
        logger.error(err)
        print(err)
        raise ValueError(err)

    if 'projectName' not in arguments:
        configuration = loadJSON("caseConfiguration.json")
        projectName = configuration['projectName']
    else:
        projectName = arguments.projectName

    logger.info(f"Using project {projectName}")

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    try:
        dispersionWorkFlow = loadJSON(arguments.dispersionWorkflow)
    except ValueError:
        dispersionWorkFlow = None

    if dispersionWorkFlow is None:
        err = f"The dispersionCaseDirectory must be a workflow file. Got {arguments.dispersionCaseDirectory}"
        logger.error(err)
        raise ValueError(err)

    # 1. Check if the dispersionFlowField is a workflow name or a directory
    dispersionFlowFieldDocumentList = tk.getWorkflowDocumentFromDB(arguments.dispersionFlowField)
    if len(dispersionFlowFieldDocumentList) == 0:
        logger.debug(f"Not found in the DB, trying to use the argument {arguments.dispersionFlowField} as a directory")
        if os.path.isdir(os.path.abspath(arguments.dispersionFlowField)):
            dispersionFlowFieldName = os.path.abspath(arguments.dispersionFlowField)
            dispersionFlowFieldDirectory = dispersionFlowFieldName
            logger.info(f"Using the dispersion flow {dispersionFlowFieldName} as name and directory. Parameters will not be available since it is not in DB")
        else:
            err = f"Cannot find the dispersion flow field {arguments.dispersionFlowField}"
            logger.error(err)
            raise ValueError(err)
    else:
        dispersionFlowFieldName = dispersionFlowFieldDocumentList[0].desc['workflowName']
        dispersionFlowFieldDirectory = dispersionFlowFieldDocumentList[0].resource
        logger.debug(f"Found the dispersion flow {dispersionFlowFieldName} in the DB. Using direcotry {dispersionFlowFieldDirectory}")

    currentDispersionName = os.path.basename(arguments.dispersionWorkflow).split(".")[0]
    logger.info(f"The dispersion workflow {currentDispersionName} and the dispersionFlowField Name {dispersionFlowFieldName} (direccotry {dispersionFlowFieldDirectory})")

    updateDispersionFlowField = False
    updateWorkflow = False
    needRewrite = False
    needDBAdd = False

    # 2. Check if the dispersion workflow name is in the DB.
    DBDispersionWorkflowDocument = tk.getWorkflowDocumentFromDB(currentDispersionName) #, dispersionFlowFieldName=dispersionFlowFieldName)
    if len(DBDispersionWorkflowDocument) > 0:
        logger.execution(f" the workflow {currentDispersionName} in the DB. First see if the dispersion flow field is similar")

        # 2.1 Check if the dispersion flow fields are similar
        if (DBDispersionWorkflowDocument[0]['desc']['dispersionFlowFieldName'] != dispersionFlowFieldName):
            logger.debug(f"The dispersion with the name {currentDispersionName} exists in DB and is associated to the disperesion flow {DBDispersionWorkflowDocument[0]['desc']['dispersionFlowFieldName']}.")
            if arguments.updateDB:
                logger.debug("updateDB and rewrite flags exist, so mark for removing the dispersion case and updating the record")
                updateDispersionFlowField = True
                needRewrite = True
            else:
                err = f"The dispersion with the name {currentDispersionName} exists in DB and is associated to the disperesion flow {DBDispersionWorkflowDocument[0]['desc']['dispersionFlowFieldName']}. In order to change it to dispersion {dispersionFlowFieldName}, either change the dispersion name or use the --updateDB"
                logger.error(err)
                raise ValueError(err)
        else:
            logger.debug("dispersion not found in DB, add it")
            needDBAdd = True

        logger.execution(f"Check if the disk workflow is identical to the DB workflow")
        diskDispersionWorkflow = tk.getHermesWorkflowFromJSON(dispersionWorkFlow)
        DBDispersionWorkflow   = tk.getHermesWorkflowFromJSON(DBDispersionWorkflowDocument[0]['desc']['workflow'])

        # 2.1 Check if the workflows are similar
        res = tk.compareWorkflowsObj([DBDispersionWorkflow, diskDispersionWorkflow])
        if len(res.columns) == 1:
            logger.info(f"The workflow in the DB is identical to the disk.")
            DispersionWorkflow = diskDispersionWorkflow
        else:
            logger.info(f"The workflow in the DB is different than the disk. Must specify whether to use the disk version  (--updateDB), or the DB version (--exportFromDB). ")
            if arguments.exportFromDB:
                logger.execution(f"Exporting the workflow to file {arguments.dispersionWorkflow}")
                with open(arguments.dispersionWorkflow, 'w') as JSONOut:
                    json.dump(DBDispersionWorkflowDocument[0]['desc']['workflow'],JSONOut,indent=4)

                DispersionWorkflow = DBDispersionWorkflow
            elif arguments.updateDB:
                logger.execution("Updating the DB with the workflow from disk")
                DispersionWorkflow = diskDispersionWorkflow
                updateWorkflow = True
                needRewrite = True
            else:
                err = f"The workflow in {arguments.dispersionWorkflow} file and in the DB are different. Use the --updateDB to update the DB or --exportFromDB to rewrite the file"
                logger.error(err)
                raise ValueError(err)

    # 3. Check if workflow already in DB under a different name
    logger.execution("Check if workflow already in DB under a different name")
    #    Stop unless allowDuplicate flag exist
    DBDispersionWorkflow = tk.getHermesWorkflowFromDB(dispersionWorkFlow, dispersionFlowFieldName=dispersionFlowFieldName)

    if DBDispersionWorkflow is not None:
        logger.debug(f"Found the workflow in the DB under the name {DBDispersionWorkflow[0]['desc']['workflowName']}")
        if currentDispersionName != DBDispersionWorkflow[0]['desc']['workflowName'] and not arguments.allowDuplicate:
            info = f"Current simulation name is {currentDispersionName}. The workflow already in DB under the name {DBDispersionWorkflow[0]['desc']['workflowName']}. use the --allowDuplicate to continue"
            logger.error(info)
            raise ValueError(info)

    dispersionName = arguments.dispersionWorkflow.split(".")[0]
    dispersionDirectoryName = os.path.abspath(dispersionName)
    logger.debug(f"Getting the dispersion directory name from {arguments.dispersionWorkflow}: using {dispersionDirectoryName}")

    # 4.  check if the dispersion directory exists.
    #    if it exist, remove if the --rewrite flag exists.
    logger.execution("Check if dispersion already exists on the disk. If it is remove for fresh start only in rewrite flag exists")
    if os.path.isdir(dispersionDirectoryName):
        if needRewrite:
            logger.execution(f"The directory {dispersionDirectoryName} exists and differ than DB. Rewriting is needed as disk version was selected with --updateDB flag")
        if arguments.rewrite:
            logger.info(f"rewrite flag exists: removing the directory {dispersionDirectoryName}")
            shutil.rmtree(dispersionDirectoryName)
        else:
            err = f"The dispersion directory {dispersionDirectoryName} already exists. use --rewrite to force removing and recreating"
            logger.error(err)
            raise ValueError(err)

    # 3. Create the dispersion case and link to the workflow.
    logger.info(f"Creating dispersion case {dispersionDirectoryName}  and linking to {dispersionFlowFieldDirectory}")
    tk.stochasticLagrangian.createDispersionCaseDirectory(dispersionDirectoryName,dispersionFlowDirectory=dispersionFlowFieldDirectory)

    # 4. Build the workflow and execute it.
    logger.info("Building the workflow and executing it, forcing rebuild and reexecute")
    arguments.force = True
    arguments.workflow = arguments.dispersionWorkflow

    handler_buildExecute(arguments)

    # 5. Add/update to the DB.
    logger.execution("add to DB if not exist. If exists, check what needs to be updated (dispersionFlowField or the entire code)")
    if needDBAdd:
        groupName = dispersionName.split("_")[0]
        groupID   = dispersionName.split("_")[1]
        logger.debug(f"Adding new document with group {groupName} and {groupID}")

        parameters = DispersionWorkflow.parametersJSON

        if len(DBDispersionWorkflowDocument) > 0:
            parameters['flowParameters'] = DBDispersionWorkflowDocument[0].desc['flowParameters']

        doc = tk.addSimulationsDocument(resource=dispersionFlowFieldDirectory,
                                          dataFormat=datatypes.STRING,
                                          type=workflowsTypes.OF_DISPERSION.value,
                                          desc=dict(
                                              dispersionFlowFieldName=dispersionFlowFieldName,
                                              groupName=groupName,
                                              groupID=groupID,
                                              workflowName=dispersionDirectoryName,
                                              workflowType=DispersionWorkflow.workflowType,
                                              workflow=DispersionWorkflow.json,
                                              parameters=parameters)
                                          )

    elif (updateDispersionFlowField or updateWorkflow):
        doc = dispersionFlowFieldDocumentList[0]
        if updateDispersionFlowField:
            logger.debug("Updating dispersion flow field")
            doc.desc['dispersionFlowFieldName'] = dispersionFlowFieldName

            if len(DBDispersionWorkflowDocument) > 0:
                doc.desc['parameters']['flowParameters'] = DBDispersionWorkflowDocument[0].desc['flowParameters']


        if updateWorkflow:
            logger.debug("Updating workflow")
            doc.desc['workflow'] = DispersionWorkflow.json
            doc.desc['parameters'] = DispersionWorkflow.parametersJSON

            if len(DBDispersionWorkflowDocument) > 0:
                doc.desc['parameters']['flowParameters'] = DBDispersionWorkflowDocument[0].desc['flowParameters']

        logger.debug("Save to DB")
        doc.save()

