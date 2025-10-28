import json
import logging
import os

import pandas
from hera.utils import jsonutils

from .. import toolkitHome
from ..simulations.hermesWorkflowToolkit import hermesWorkflowToolkit
from ..utils import loadJSON,compareJSONS
from hermes import workflow
from hermes.utils.workflowAssembly import handler_build,handler_execute


def WorkflowsGroup_list(args):
    """
        List all the simulation groups
    Parameters
    ----------
    args
        projectName : the project name. If not supplied, get from caseConfiguration

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin.hera_workflows.group_list")
    logger.info(" -- Starting: listing all the simulation groups --")

    _confirm_project_name(args, logger)

    logger.info(f"Using project name {args.projectName}. ")


    solver = args.solver
    workflowName = args.workflowName

    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=args.projectName)
    wftk.listGroups(solver=solver, workflowName=workflowName)


def workflow_add(args):
    """
        Adds a new simulation to the group. The workflow is expanded and store in the db under the requested group.

        If the groupName is None, extract the group name from the simulation name. that is,
        assume that the simulation name is <group name>_<id>.

    Parameters
    ----------
    args:
        workflow : The workflow file that will be prepared.
        groupName : The group of the simulation.
        variations : Json with the parameter variations on the requested workflow.
                    The structure of the JSON is:
                    {
                        "parameterVariation" : {
                                TBD
                        }
                    }
        simulationType: the type of the simulation.

        overwrite : bool,
                overwite thw workflow with the given name
        force:   bool,
            add the simulation to the db, even if the workflow exists under a different name.
        assignName: bool,
            generate automated name to the workflow

    Returns
    -------
        None
    """
    logger = logging.getLogger("hera.bin.hera_workflows.add")
    logger.info(" -- Starting: adding workflow to the group --")

    _confirm_project_name(args, logger)
    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=args.projectName)

    workflowFile = args.workflow
    wftk.addUpdateWorkflowFileInGroup(workflowFile)


def workflow_delete(arguments):
    """
        deletes the simulation/list of simulations.old from the DB.
        The default is to export them to the disk (unless no-export flag is supplied
    Parameters
    ----------
    arguments:
        projectName: if not supplied get from the deleted object.

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin.hera_workflows.delete")
    logger.info(f" -- Starting: Deleting workflows --")

    _confirm_project_name(arguments, logger)
    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=arguments.projectName)

    simulationList = wftk.getWorkflowListDocumentFromDB(list(arguments.workflows))
    completeRemove = []
    for sim in simulationList:
        shouldRemove = True
        logger.info(f" Deleting the workflow: {sim['desc']['workflowName']}")
        outfileName = f"{sim['desc']['workflowName']}.json"

        if arguments.Export:
            logger.debug(f"Exporting the deleted document as {outfileName}")
            if not os.path.isfile(outfileName) or arguments.forceOverwrite:
                with open(outfileName,"w") as outfile:
                    outjson = dict(workflow=sim['desc']['workflow'])
                    json.dump(outjson,outfile,indent=4)
            else:
                logger.debug(f"...workflow {sim['desc']['workflowName']} (file {outfileName}) exists in current directory. Skipping Remove. To enforce removing either use the no-export or the forceOverwrite flags")
                shouldRemove = False

        if shouldRemove:
            res = sim.resource
            completeRemove.append(f"shutil.rmtree('{res}')")
            logger.debug("... remove from DB")
            sim.delete()

    with open("completeDelete.py","w") as outfile:
        outfile.write("\n".join(completeRemove))

    print("In order to remove all directories from disk type: python completeRemove.py")

def workflow_export(arguments):
    """
        Exports the workflow in the DB to the disk

    Parameters
    ----------
    arguments:
        projectName: if not supplied get from the deleted object.

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin.hera_workflows.export")
    logger.info(f" -- Starting: Deleting workflows --")


    _confirm_project_name(arguments, logger)
    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=arguments.projectName)
    simulationList = wftk.getWorkflowListDocumentFromDB(list(arguments.workflows))

    for sim in simulationList:
        outfileName = f"{sim['desc']['workflowName']}.json"

        logger.debug(f"Exporting document as {outfileName}")
        if not os.path.isfile(outfileName) or arguments.forceOverwrite:
            with open(outfileName,"w") as outfile:
                json.dump(sim['desc']['workflow'],outfile,indent=4)
        else:
            logger.debug(f"...workflow {sim['desc']['workflowName']} (file {outfileName}) exists in current directory, not export. Removing file or use the forceOverwrite flags")
            print(f"...workflow {sim['desc']['workflowName']} (file {outfileName}) exists in current directory, not export. Removing file or use the forceOverwrite flags")

def workflow_compareToDisk(arguments):
    """
        Exports the workflow in the DB to the disk

    Parameters
    ----------
    arguments:
        projectName: if not supplied get from the deleted object.

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin.hera_workflows.compareToDisk")
    logger.info(f" -- Starting: Deleting workflows --")


    _confirm_project_name(arguments, logger)
    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=arguments.projectName)
    simulationList = wftk.getHermesWorkflowFromDB(list(arguments.workflows),returnFirst=False)

    for sim in simulationList:
        outfileName = f"{sim.name}.json"
        if not os.path.isfile(outfileName):
            print(f"Workflow {sim.name} (file {outfileName} does not exist on the disk. use export to create it. ")
        else:
            localWorkflow = wftk.getHermesWorkflowFromJSON(loadJSON(outfileName),name="Local")
            smName = sim.name
            sim.name = "DB"
            res = compareJSONS(DB=sim.parametersJSON,LocalFile=localWorkflow.parametersJSON)
            ttl = f"Simulation {smName}"
            print(ttl)
            print("-"*len(ttl))
            if res.empty:
                print("\t\t ** Disk and DB are identical")
                print(" ")
            else:
                print(res)

def sorround_with_char(text:str, total_len:int, char:str="-"):
    chars = char*max(((total_len-len(text))//2 ), 0)
    return chars + text + chars

def workflow_sync_to_db(arguments):
    """Check for differences in the local workflow vs DB and update the changes

    Parameters
    ----------
    arguments:
        projectName: if not supplied get from the deleted object.

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin.hera_workflows.sync_to_db")
    logger.info(f" -- Starting: syncing workflows --")

    _confirm_project_name(arguments, logger)

    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=arguments.projectName)
    assert isinstance(wftk, hermesWorkflowToolkit)
    workflowDocList = wftk.getWorkflowListDocumentFromDB(list(arguments.workflows))

    if len(workflowDocList) == 0:
        logger.error(f"... not found.")

    for workflowDoc in workflowDocList:
        sim = wftk.getHemresWorkflowFromDocument(documentList=workflowDoc)
        outfileName = f"{sim.name}.json"
        if not os.path.isfile(outfileName):
            print(f"Workflow {sim.name} (file {outfileName}) does not exist on the disk. use export to create it.")
        else:
            localWorkflow = wftk.getHermesWorkflowFromJSON(loadJSON(outfileName),name="Local")
            simName = sim.name
            sim.name = "DB"
            res = compareJSONS(DB=sim.parametersJSON,LocalFile=localWorkflow.parametersJSON)
            assert isinstance(res, pandas.DataFrame)
            title = sorround_with_char(f"Simulation {simName}", 64)
            if not arguments.quiet:
                print(title)
            
            if not arguments.quiet:
                print("Disk and DB parameters are identical" if res.empty else f"Found Changes:\n{res}")
            if (not res.empty) or arguments.force:
                logger.info(f"Updating DB with the changes for {sim.name}")
                wftk.updateDocumentWorkflow(document=workflowDoc, json=outfileName)

def create_workflow_variations(arguments):
    """creates variations to the base workflow provided based on the variation json

    Parameters
    ----------
    arguments:
        projectName: if not supplied get from the deleted object.

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin.hera_workflows.create_workflow_variations")
    logger.info(f" -- Starting: creating workflow variation --")
    _confirm_project_name(arguments, logger)
    variation_json_path = str(arguments.variationJson)
    project_name = str(arguments.projectName)
    workflow_name = str(arguments.workflow)
    dry_run = bool(arguments.dry_run)
    overwrite = bool(arguments.overwrite)
    naming_scheme = str(arguments.naming_scheme)


    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=project_name)
    assert isinstance(wftk, hermesWorkflowToolkit)
    workflow_doc_list = wftk.getWorkflowListDocumentFromDB(workflow_name)

    if len(workflow_doc_list) == 0:
        logger.error("workflow not found.")
    wokflow_doc = workflow_doc_list[0]
    base_workflow_json = wokflow_doc.to_mongo().copy()['desc']['workflow']
    base_workflow_name = str(wokflow_doc['desc']['groupName'])
    base_workflow_path = str(wokflow_doc['resource'])
    
    if not os.path.exists(variation_json_path):
        logger.error(f"{variation_json_path} doesn't point to anything.")    
    elif not os.path.isfile(variation_json_path):
        logger.error(f"{variation_json_path} must be a file")
    
    with open(variation_json_path) as varitation_json_file:
        variation_json = json.load(varitation_json_file)
    
    variation_scheme = None
    # Might be changed to a dictionary {<scheme>: <scheme_func>}
    if naming_scheme == "index":
        variation_scheme = enumerate(jsonutils.JSONVariations(base_workflow_json, variation_json))
    else:
        logger.error(f"naming scheme {naming_scheme} is not implemented")
        return
    
    for variation_suffix, variation in variation_scheme:
        variation_path = os.path.join(base_workflow_path, "..", base_workflow_name+"_"+str(variation_suffix))
        if os.path.exists(variation_path) and not overwrite:
            logger.warning(f"skipping variation {variation_suffix} since it already exists")
        else:
            with open(variation_path, "w") as variation_file:
                json.dump(variation, variation_file)
        if not dry_run:
            _ = wftk.addUpdateWorkflowFileInGroup(variation_path)

def batch_delete_workflows(arguments):
    """creates variations to the base workflow provided based on the variation json

    Parameters
    ----------
    arguments:
        projectName: if not supplied get from the deleted object.

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin.hera_workflows.batch_delete_workflows")
    logger.info(f" -- Starting: batch delete workflows --")
    _confirm_project_name(arguments, logger)
    project_name = str(arguments.projectName)
    workflow_name = str(arguments.workflow) if arguments.workflow else None
    group_name = str(arguments.groupName) if arguments.groupName else None
    hard_delete = bool(arguments.hard)

    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=project_name)
    assert isinstance(wftk, hermesWorkflowToolkit)
    
    workflow_to_exclude = None
    if workflow_name:
        workflow_doc_list = wftk.getWorkflowListDocumentFromDB(workflow_name)

        if len(workflow_doc_list) == 0:
            raise ValueError("workflow not found.")
        else:
            workflow_to_exclude=workflow_doc_list[0]
    elif group_name is None:
        raise ValueError("either group name or workflow name must be provided")

    group_name = group_name if group_name else workflow_to_exclude['desc']['groupName']
    wftk.deleteWorkflowInGroup(group_name,hard_delete,resetCounter=True, exclude=[workflow_name])


def _confirm_project_name(arguments, logger):
    if arguments.projectName is None:
        logger.debug(
            f"projectName is not provided. Looking for the project name in the caseConfiguration.json file (projectName key) ")
        caseConfiguration = loadJSON("caseConfiguration.json")
        arguments.projectName = caseConfiguration['projectName']



def workflow_list(arguments):
    """
            Lists the simulations of the project.

            Parameters
            ----------
            arguments : argument struct with the field:

                - object : str, optional
                            List only the group name if exists.
.
                - projectName: str, optional
                            If not supplied, use the project name in caseConfiguration.json.

                - workflowGroup: str, optional
                            If not supplied, use the file name.

                - no-nodes    : bool [default : true].
                            If exists, does not list the node names.
                            if --parameters exists, then this option is ignored.

                - parameters [node list]: list of string, optional
                            If flag exists lists only the parameters that are different between different simulations.
                            unless --all exists. In that case, list all nodes.

                            if [node list] is not empty, lists all the parameters of this node.

                - all : boolean, [optional]
                        If --parameters flag is used, lists all the parameters if exists. Otherwise
                        list only the nodes that were changed.

            Returns
            -------
                    A string.

                    * without --parameters flag:

                    <work-flow type>
                    ----------------
                            + <group name>
                                    - simulation name : [node list]
                                    .
                                    .


                    * with --parameters exists:

                    <work-flow type> - <group name>
                     |   simulation name |  <parameters>   |
                     +-------------------+-----------------+
    """
    logger = logging.getLogger("hera.bin.hera_workflows.list")
    logger.info(f" -- Starting: Listing simulations --")

    _confirm_project_name(arguments, logger)

    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=arguments.projectName)

    # if arguments.object is None:
    #     # Listing all the groups in the toolkit.
    #     docList =  wftk.getSimulationsDocuments(type=wftk.WORKFLOW)
    #     if docList is None:
    #         logger.info(f"There are no hermes workflows in project {projectName}")
    #
    #     groupNameList = set([x['desc']['groupName'] for x in docList])
    #
    #     title = f"The simulation groups in project *{projectName}* "
    #     print(title)
    #     print("-" * len(title))
    #     print("\n".join([x for x in groupNameList]))
    #
    # else:

    simDocument = wftk.getWorkflowListDocumentFromDB(arguments.group)
    if len(simDocument) == 0:
        print(f"{arguments.object} is not a simulation, directory, workflow file or a simulation group in project {arguments.projectName} ")
    workflowGroup = simDocument[0].desc[wftk.DESC_GROUPNAME]

    listNodes     = arguments.nodes
    parameters    = arguments.parameters

    simulationList = wftk.listWorkflows(workflowGroup=workflowGroup,
                                        listNodes=listNodes,
                                        listParameters=parameters)


    title = f"The simulations in group *{workflowGroup}*  in project *{arguments.projectName}* "
    print(title)
    print("-" * len(title))

    for doc in simulationList:
        print(f"\t* {doc['workflowName']}")
        if listNodes:
                for node in doc['nodes']:
                    print(f"\t\t + {node}")

        if parameters:
            for nodeName,nodeData in doc['parameters'].items():

                print(f"\t\t + {nodeName}")
                for pname,pvalue in nodeData.items():
                    print(f"\t\t\t  - {pname}")


def workflowNodes_list(arguments):
    """
        Lists the nodes in the requested workflow. The workflow can be a file on the disk or a name
        of a simulation in the database.

    Parameters
    ----------
    arguments
            projectName: str
                The name of the project.
            workflowName: str
                A file on the disk or a simulation in the DB.
    Returns
    -------
        prints a list of all the nodes of the workflow.
    """
    logger = logging.getLogger("hera.bin.hera_workflows.listNodes")
    logger.info(f" -- Starting: Listing workflow nodes --")

    if os.path.exists(arguments.workflowName) and not os.path.isdir(arguments.workflowName):
            json = loadJSON(arguments.workflowName)
            hermesObject = workflow(json)
    else:
        if arguments.projectName is None:
            raise ValueError("Must supply a project name for a non-file workflow")

        wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=arguments.projectName)

        hermesObject = wftk.getHermesWorkflowFromDB(arguments.workflowName)

    tlte = f"The nodes of the {arguments.workflowName}"
    print(tlte)
    print("-"*len(tlte))

    if arguments.parameters:
        for hnodeName,hnodeData in hermesObject.items():
            print(f"\t * {hnodeName}")
            for prop in hnodeData.parameters.keys():
                print(f"\t\t + {prop}")
    else:
        print("\t * "+"\n\t * ".join(hermesObject.nodeList))

def workflowNodes_listParameters(arguments):
    """
            List the parameters of the node.

    Parameters
    ----------
    arguments
            projectName: str
                The name of the project.
            nodename : str
                The name of the node to list.

            workflowName: str
                A file on the disk or a simulation in the DB.

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin.hera_workflows.listNodeParameters")
    logger.info(f" -- Starting: Listing node parameters --")

    if os.path.isfile(arguments.workflowName):
            json = loadJSON(arguments.workflowName)
            hermesObject = workflow(json)
    else:
        wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=arguments.projectName)

        hermesObject = wftk.getHermesWorkflowFromDB(arguments.workflowName)

    if arguments.nodeName not in hermesObject.nodeList:
        raise ValueError(f" Node {arguments.nodeName} not found in workflow {arguments.workflowName}. Existing nodes are: {','.join(hermesObject.nodeList)}")

    tlte = f"The parameters of the node {arguments.nodeName} in the workflow {arguments.workflowName}"
    print(tlte)
    print("-"*len(tlte))
    import json

    for nd,pm in hermesObject[arguments.nodeName].parameters.items():
        vls = json.dumps(pm,indent=4)
        print(f"-\t {nd}:  {vls}")

def workflow_compare(arguments):
    """
            Compares the parameters of the list of simulations that were supplied.
            Specifically for a list of simulations. Comparing all the simulations of a group is
            achieved with list simulations.

    Parameters
    ----------
    arguments
        projectName : str, the name of the projet

        simulations: [groupName] - compare all the simulations.old,
                     [sim1,sim2,..] compare the different simulations.old. simX is either a simulation name in the DB or a file on the disk

        format : The format of the output.
            - pandas
            - latex
            - csv
            - json

        file : None
            If not None, save the output to the file.

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin.hera_workflows.compare")
    logger.info(f" -- Starting: Listing workflow nodes --")

    _confirm_project_name(arguments, logger)
    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=arguments.projectName)

    res = wftk.compareWorkflows(arguments.workflows, longFormat=arguments.longFormat, transpose=arguments.transpose)


    if arguments.format == "pandas":
        output = res
        ext = "txt"
    elif arguments.format == "latex":
        output = res.to_latex()
        ext = "tex"
    elif arguments.format == "csv":
        output = res.to_csv()
        ext = "csv"
    else:
        output = json.dumps(loadJSON(res.to_json()),indent=4)
        ext = "json"

    if len(res)==0:
        print(f"Could not found any workflows to compare in project {arguments.projectName}")
    else:
        print(output)

        if arguments.file is not None:
            flName = arguments.file if "." in arguments.file else f"{arguments.file}.{ext}"

            with open(flName,"w") as outputFile:
                outputFile.write(output)

def workflow_build(arguments):
    handler_build(arguments)

def workflow_execute(arguments):
    handler_execute(arguments)

def workflow_buildExecute(arguments):
    """
        Builds the executes the workflow.

    Parameters
    ----------
    arguments
        workflowName: str
            The name of the workflow in the DB or the file name.

        addDB : bool
            If true, add the workflow if not found.

    Returns
    -------

    """
    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=arguments.projectName)

    workflowName = arguments.workflowName.split(".")[0]

    docList = wftk.getWorkflowListDocumentFromDB(workflowName)
    if len(docList) == 0:
        if os.path.isfile(arguments.workflowName):
            wftk.addUpdateWorkflowFileInGroup(arguments.workflowName)
        else:
            raise ValueError(f"The workflowName {arguments.workflowName} is not in the DB and not a file on the disk")

    wftk.executeWorkflowFromDB(workflowName)