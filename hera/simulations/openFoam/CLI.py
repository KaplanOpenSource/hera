import json
import logging
import glob
import os
import shutil
import numpy
import pandas

from ...datalayer import datatypes
from ... import toolkitHome
from ...utils.jsonutils import loadJSON
from ...utils.freeCAD import getObjFileBoundaries
from ...utils.logging import get_logger


def simpleFoam_createEmpty(arguments):
    logger = logging.getLogger("hera.bin")
    logger.execution(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    if 'projectName' not in arguments:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"
        try:
            configuration = loadJSON(configurationFile)
        except:
            err = f"Configuration file {configurationFile} not found! creating a basic file"
            with open(configurationFile,'w') as defaultConfFile:
                json.dump(dict(projectName=None),defaultConfFile,indent=4)

            raise ValueError(err)

        projectName = configuration['projectName']
    else:
        projectName = arguments.projectName

    logger.info(f"Using project {projectName}")
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    
    tk.createEmptyCase(caseDirectory = arguments.caseDirectory,
                       fieldList = arguments.fields,
                       simulationType=tk.SIMULATIONTYPE_INCOMPRESSIBLE,
                       additionalFieldsDescription = arguments.fieldsDescription)


    logger.execution(f"----- End -----")

def simpleFoam_templates_list(arguments):
    logger = logging.getLogger("hera.bin")
    logger.execution(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    projectName = None if 'projectName' not in arguments else arguments.projectName # from hera 2.13.2 the toolkit searches the project name in the case file.
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    templates = tk.listHermesTemplates(arguments.solver)

    ttl = f"The templates for project {tk.projectName} with solver {arguments.solver}"
    print()
    print("-" * len(ttl))
    print(ttl)
    print("-"*len(ttl))
    print(templates)

def simpleFoam_templates_create(arguments):
    logger = logging.getLogger("hera.bin")
    logger.execution(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=None)

    projectName = arguments.projectName # from hera 2.13.2 the toolkit searches the project name in the case file.
    outputPath = projectName if arguments.projectPath is None else arguments.projectPath
    tk.createProjectDirectory(outputPath=outputPath,projectName=projectName)

    templates = tk.listHermesTemplates(arguments.solver)
    if arguments.templateName not in templates.index:
        err = f"{arguments.templateName} is not known. Use one of the \n " + str(templates)
        print(err)
        raise ValueError(err)

    groupName = arguments.templateName if arguments.groupName is None else arguments.groupName

    with open(os.path.join(outputPath,f"{groupName}_1.json"),"w") as outFile:
        json.dump(tk.getDatasourceData(arguments.templateName),outFile, indent=4)


def stochasticLagrangian_dispersionFlow_create(arguments):
    logger = logging.getLogger("hera.bin")
    logger.execution(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    if ('projectName' in arguments) and (arguments.projectName is not None):
        projectName = arguments.projectName
    else:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"
        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']

    logger.info(f"Adding dispersion flow to project {projectName}")
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)

    if ('DFF' in arguments) and (arguments.DFF is not None):
        dffList = arguments.DFF
    else:
        dffList = [x for x in configuration["DispersionFlows"].keys()]

    logger.info(f"Createing dispersion flows {','.join(dffList)}")

    for flowName in dffList:
        logger.debug(f"Creating dispersion flow : {flowName}")
        try:
            flowdata = configuration["DispersionFlows"][flowName]
            tk.stochasticLagrangian.createDispersionFlowField(flowName=flowName,flowData=flowdata,OriginalFlowField=arguments.OriginalFlowField,overwrite=arguments.overwrite)
        except KeyError:
            err = f"Flow field {flowName} not Found. Found the following flows: {','.join(configuration['flowFields']['Flows'].keys())}"
            logger.error(err)
            raise KeyError(err)
        except FileExistsError:
            err = f"Flow field {flowName} Already exists. Use --overwrite to recreate"
            logger.error(err)

    logger.execution(f"----- End -----")


def stochasticLagrangian_dispersionFlow_list(arguments):
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
    res = tk.workflowCompare(workflowsType="stochasticLagrangianSolver")

    if arguments.format == "pandas":
        for group in res.groupName.unique():
            ttl = f"Group name {group}"
            print(ttl)
            print("-" * len(ttl))
            print(res.query('groupName == @group'))

    elif arguments.format == "latex":
        output = res.to_latex()
        ext = "tex"
    elif arguments.format == "csv":
        output = res.to_csv()
        ext = "csv"
    else:
        output = json.dumps(loadJSON(res.to_json()), indent=4)
        ext = "json"

    if len(res) == 0:
        print(f"Could not found any workflows to compare in project {projectName}")
    else:
        print(output)

        if arguments.file is not None:
            flName = arguments.file if "." in arguments.file else f"{arguments.file}.{ext}"

            with open(flName, "w") as outputFile:
                outputFile.write(output)


def stochasticLagrangian_dispersion_create(arguments):
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
    logger = logging.getLogger("hera.bin.stochasticLagrangian.dispersion_create")
    logger.execution(f"----- Start -----")
    logger.debug(f"  Got arguments: {arguments}")

    if ('projectName' in arguments) and (arguments.projectName is not None):
        projectName = arguments.projectName
    else:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"
        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']

    logger.info(f"Using project {projectName}")

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)

    dispersionDirectoryName = arguments.dispersionName
    logger.info(f"Creating dispersion case {dispersionDirectoryName}")

    dispersionFlowFieldName = arguments.dispersionFlowField
    logger.info(f"Getting the dispersion flowField {dispersionFlowFieldName} ")
    doc = tk.getWorkflowDocumentFromDB(dispersionFlowFieldName, tk.OF_FLOWDISPERSION)
    if len(doc)==0:
        logger.info(f"Dispersion flow {dispersionFlowFieldName} not found in DB. Trying to use as a directory")
        if not os.path.exists(dispersionFlowFieldName):
            err = f"{dispersionFlowFieldName} not found!. "
            logger.error(err)
            raise ValueError(err)
        dispersionFlowField = os.path.abspath(dispersionFlowFieldName)
    else:
        logger.info("Found the dispersion in the DB. Using it")
        dispersionFlowField = doc[0].resource

    if os.path.exists(dispersionDirectoryName):
        logger.info(f"Dispersion {dispersionDirectoryName} exists.")
        if arguments.overwrite:
            logger.info("Got the overwrite flag: removing old directory")
            shutil.rmtree(dispersionDirectoryName,)
        else:
            err = "Dispersion flow exists. Use --overwrite to force recreation"
            raise ValueError(err)

    # 3. Create the dispersion case and link to the workflow.
    logger.info(f"Creating dispersion case {dispersionDirectoryName}  and linking to {dispersionFlowField}")
    tk.stochasticLagrangian.createAndLinkDispersionCaseDirectory(dispersionDirectoryName,dispersionFlowDirectory=dispersionFlowField)


def stochasticLagrangian_source_cylinder(arguments):
    logger = logging.getLogger("hera.bin.stochasticLagrangian_source_cylinder")
    logger.execution(f"----- Start -----")
    logger.debug(f"  Got arguments: {arguments}")

    center = arguments.center
    params = dict(x=float(center[0]),
                  y=float(center[1]),
                  z=float(center[2]),
                  radius = float(arguments.radius),
                  height = float(arguments.height),
                  nParticles=int(arguments.particles)
                  )
    if ('projectName' in arguments) and (arguments.projectName is not None):
        projectName = arguments.projectName
    else:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"
        logger.debug(f"Loading configuration file {configurationFile}")
        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']

    dispersionName = arguments.dispersionName

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    tk.stochasticLagrangian.writeParticlePositionFile(type="Cylinder",dispersionName=dispersionName, **params)


def stochasticLagrangian_source_makeEscapedMassFile(args):
    """
        Creates an openfoam mass overtime file for the lagrangian output
        from a dataframe.

    Parameters
    ----------
    args

    Returns
    -------

    """

    case   = os.path.abspath(args.casePath)
    massFileName = f"{args.patch}Mass" if args.massFileName is None else args.massFileName
    dt = args.dt
    LSMtoolkit = toolkitHome.getToolkit(toolkitHome.OF_LSM,"tmpProject",casePath=case)
    data = LSMtoolkit.analysis.getMassFromLog(logFile=args.logFile,solver=args.solver)
    data = data.loc[data.filterType == args.patch].loc[data.action == args.action]
    data["diffMass"] = data.mass.diff()
    data = data.fillna(0)
    if dt is None:
        timesteps = data.time
    else:
        timesteps = numpy.arange(data.time.min(),data.time.max(),float(dt))
        times = pandas.DataFrame({"time":timesteps})
        data = data.set_index("time").join(times.set_index("time"),how="outer").reset_index()
        data = data.interpolate()
    newstr = "/*--------------------------------*- C++ -*----------------------------------*\n" \
             "| =========                 |                                                 |\n" \
             "| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n" \
             "|  \    /   O peration     | Version:  dev                                   |\n" \
             "|   \  /    A nd           | Web:      www.OpenFOAM.org                      |\n" \
             "|    \/     M anipulation  |                                                 |\n" \
             "\*---------------------------------------------------------------------------*/\n" \
             "FoamFile\n{    version     2.0;\n    format      ascii;\n    class       scalarField;\n    object      kinematicCloudPositions;\n}\n" \
             f"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n{len(data)}\n(\n"
    for time in timesteps:
        newstr += f"{float(data.loc[data.time==time].diffMass)}\n"
    newstr += ")"
    print("saving in ",os.path.join(case,"constant",massFileName))
    f = open(os.path.join(case,"constant",massFileName), "w")
    f.write(newstr)
    f.close()

def stochasticLagrangian_postProcess_toParquet(arguments):
    """
        Creates a parquet file from the case results
    """
    logger = logging.getLogger("hera.bin.stochasticLagrangian_postProcess_toParquet")
    logger.execution(f"----- Start -----")
    logger.debug(f"  Got arguments: {arguments}")

    if ('projectName' in arguments) and (arguments.projectName is not None):
        projectName = arguments.projectName
    else:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"
        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']

    logger.info(f"Using project {projectName} for cloud name {arguments.cloudName}")

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    data = tk.stochasticLagrangian.getCaseResults(caseDescriptor=arguments.case,withVelocity=True,withMass=True,overwrite=arguments.overwrite,cloudName=arguments.cloudName)
    tk.stochasticLagrangian.presentation.toVTU(data,outDirectory=arguments.outputDirectory,outFile=os.path.basename(arguments.case))

def stochasticLagrangian_postProcess_toVTK(arguments):
    """
        Creates VTK files from the case.
        Creates the parquet file if does not exist.
    """
    logger = logging.getLogger("hera.bin.stochasticLagrangian_postProcess_toParquet")
    logger.execution(f"----- Start -----")
    logger.debug(f"  Got arguments: {arguments}")

    if ('projectName' in arguments) and (arguments.projectName is not None):
        projectName = arguments.projectName
    else:
        configurationFile = arguments.configurationFile if 'configurationFile'  in arguments else "caseConfiguration.json"
        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']

    logger.info(f"Using project {projectName}")

    if 'outputdir' in arguments:
        base_outputdir = os.getcwd() if arguments.outputdir is None else arguments.outputdir
    else:
        base_outputdir = os.getcwd()

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    docList = tk.getWorkflowDocumentFromDB(arguments.case)
    if len(docList) == 0:
        logger.info(f"The name {arguments.case} is not found. Try to use it as a directory")
        if os.path.isdir(arguments.case):
            logger.info(f"Found {arguments.case} as directory. Trying to use it as lagrangian simulation and save it in VTK format")
            outputName = arguments.case
            cache = False
    else:
        doc = docList[0]
        outputName = doc.desc['workflowName']
        cache = True

    outputdir = os.path.join(base_outputdir,"VTK",outputName)
    logger.info(f"Writing values to directory {outputdir}")
    os.makedirs(outputdir,exist_ok=True)
    data = tk.stochasticLagrangian.getCaseResults(caseDescriptor=outputName,withVelocity=True,withMass=True,cache=cache)
    tk.presentation.toUnstructuredVTK(data=data, outputdirectory=outputdir, filename=arguments.cloudName, timeNameOutput=True, xcoord="x", ycoord="y", zcoord="z")


def objects_createVerticesAndBoundary(arguments):
    logger = logging.getLogger("hera.bin")
    logger.execution(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")
    try:
        import FreeCAD
        import Mesh
    except ImportError:
        logger.error("freecad module  is not installed in the environment")
        raise ImportError("freecad is not installed. please install it before trying again.")

    # Load the file
    fileName = arguments.objectFile
    Mesh.open(fileName)
    objFile  = FreeCAD.getDocument("Unnamed")

    ret = dict()
    for fieldName in arguments.fields:
        boundary = dict()
        for  regionObj in objFile.findObjects():
            boundary[regionObj.Name] = dict(type="zeroGradient")
        ret[fieldName] = dict(boundaryField=boundary)

    print("----- Bounding box vertices -----------")
    corners = getObjFileBoundaries(fileName)

    xList = ['XMin','XMax','XMax','XMin','XMin','XMax','XMax','XMin']
    yList = ['YMin','YMin','YMax','YMax','YMin','YMin','YMax','YMax']
    zList = ['ZMin','ZMin','ZMin','ZMin','ZMax','ZMax','ZMax','ZMax']

    verticesList = []
    for x,y,z in zip(xList,yList,zList):
        verticesList.append([corners[x], corners[y], corners[z]])
    print(json.dumps(dict(vertices=verticesList), indent=4))

    print("\n\n\n")
    print("----- Boundary conditions -------------")
    print(json.dumps(ret, indent=4, sort_keys=True))




# def stochasticLagrangian_dispersion_create(arguments):
#     """
#         Prepares the dispersion case:
#
#         1. Create the dispersion case (and link to dispersionFlow)
#         2. build and execute the dispersion workflow:
#             - If the workflow on disk  mismatch the DB stop unless
#                 the --exportDB or the --updateDB flags exist.
#             - If --exportDB flag exists, use the workflow in the DB.
#         3. Add to the DB:
#             - If the workflow mismatch the DB, update the DB if the --updateDB flag exists.
#             - add the parameters of the flowcase to the DB.
#
#
#
#     :param arguments:
#              - dispersionCaseDirectory : the new case of the directory
#              - dispersionFlow : the new case of the directory.
#                                 can be either a directory or a name of a dispersion flow field from the DB
#
#     :return:
#     """
#     logger = logging.getLogger("hera.bin.stochasticLagrangian.createDispersionFlow")
#     logger.execution(f"----- Start -----")
#     logger.debug(f"  Got arguments: {arguments}")
#
#     if (arguments.updateDB and arguments.exportFromDB):
#         err = "Cannot use both --updateDB and --exportFromDB"
#         logger.error(err)
#         print(err)
#         raise ValueError(err)
#
#     if 'projectName' not in arguments:
#         configuration = loadJSON("caseConfiguration.json")
#         projectName = configuration['projectName']
#     else:
#         projectName = arguments.projectName
#
#     logger.info(f"Using project {projectName}")
#
#     tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
#     try:
#         dispersionWorkFlow = loadJSON(arguments.dispersionWorkflow)
#     except ValueError:
#         dispersionWorkFlow = None
#
#     if dispersionWorkFlow is None:
#         err = f"The dispersionCaseDirectory must be a workflow file. Got {arguments.dispersionCaseDirectory}"
#         logger.error(err)
#         raise ValueError(err)
#
#     # 1. Check if the dispersionFlowField is a workflow name or a directory
#     dispersionFlowFieldDocumentList = tk.getWorkflowDocumentFromDB(arguments.dispersionFlowField)
#     if len(dispersionFlowFieldDocumentList) == 0:
#         logger.debug(f"Not found in the DB, trying to use the argument {arguments.dispersionFlowField} as a directory")
#         if os.path.isdir(os.path.abspath(arguments.dispersionFlowField)):
#             dispersionFlowFieldName = os.path.abspath(arguments.dispersionFlowField)
#             dispersionFlowFieldDirectory = dispersionFlowFieldName
#             logger.info(f"Using the dispersion flow {dispersionFlowFieldName} as name and directory. Parameters will not be available since it is not in DB")
#         else:
#             err = f"Cannot find the dispersion flow field {arguments.dispersionFlowField}"
#             logger.error(err)
#             raise ValueError(err)
#     else:
#         dispersionFlowFieldName = dispersionFlowFieldDocumentList[0].desc['workflowName']
#         dispersionFlowFieldDirectory = dispersionFlowFieldDocumentList[0].resource
#         logger.debug(f"Found the dispersion flow {dispersionFlowFieldName} in the DB. Using direcotry {dispersionFlowFieldDirectory}")
#
#     currentDispersionName = os.path.basename(arguments.dispersionWorkflow).split(".")[0]
#     logger.info(f"The dispersion workflow {currentDispersionName} and the dispersionFlowField Name {dispersionFlowFieldName} (direccotry {dispersionFlowFieldDirectory})")
#
#     updateDispersionFlowField = False
#     updateWorkflow = False
#     needRewrite = False
#     needDBAdd = False
#
#     # 2. Check if the dispersion workflow name is in the DB.
#     DBDispersionWorkflowDocument = tk.getWorkflowDocumentFromDB(currentDispersionName) #, dispersionFlowFieldName=dispersionFlowFieldName)
#     if len(DBDispersionWorkflowDocument) > 0:
#         logger.execution(f" the workflow {currentDispersionName} in the DB. First see if the dispersion flow field is similar")
#
#         # 2.1 Check if the dispersion flow fields are similar
#         if (DBDispersionWorkflowDocument[0]['desc']['dispersionFlowFieldName'] != dispersionFlowFieldName):
#             logger.debug(f"The dispersion with the name {currentDispersionName} exists in DB and is associated to the disperesion flow {DBDispersionWorkflowDocument[0]['desc']['dispersionFlowFieldName']}.")
#             if arguments.updateDB:
#                 logger.debug("updateDB and rewrite flags exist, so mark for removing the dispersion case and updating the record")
#                 updateDispersionFlowField = True
#                 needRewrite = True
#             else:
#                 err = f"The dispersion with the name {currentDispersionName} exists in DB and is associated to the disperesion flow {DBDispersionWorkflowDocument[0]['desc']['dispersionFlowFieldName']}. In order to change it to dispersion {dispersionFlowFieldName}, either change the dispersion name or use the --updateDB"
#                 logger.error(err)
#                 raise ValueError(err)
#         else:
#             logger.debug("dispersion not found in DB, add it")
#             needDBAdd = True
#
#         logger.execution(f"Check if the disk workflow is identical to the DB workflow")
#         diskDispersionWorkflow = tk.getHermesWorkflowFromJSON(dispersionWorkFlow)
#         DBDispersionWorkflow   = tk.getHermesWorkflowFromJSON(DBDispersionWorkflowDocument[0]['desc']['workflow'])
#
#         # 2.1 Check if the workflows are similar
#         res = tk.compareWorkflowsObj([DBDispersionWorkflow, diskDispersionWorkflow])
#         if len(res.columns) == 1:
#             logger.info(f"The workflow in the DB is identical to the disk.")
#             DispersionWorkflow = diskDispersionWorkflow
#         else:
#             logger.info(f"The workflow in the DB is different than the disk. Must specify whether to use the disk version  (--updateDB), or the DB version (--exportFromDB). ")
#             if arguments.exportFromDB:
#                 logger.execution(f"Exporting the workflow to file {arguments.dispersionWorkflow}")
#                 with open(arguments.dispersionWorkflow, 'w') as JSONOut:
#                     json.dump(DBDispersionWorkflowDocument[0]['desc']['workflow'],JSONOut,indent=4)
#
#                 DispersionWorkflow = DBDispersionWorkflow
#             elif arguments.updateDB:
#                 logger.execution("Updating the DB with the workflow from disk")
#                 DispersionWorkflow = diskDispersionWorkflow
#                 updateWorkflow = True
#                 needRewrite = True
#             else:
#                 err = f"The workflow in {arguments.dispersionWorkflow} file and in the DB are different. Use the --updateDB to update the DB or --exportFromDB to rewrite the file"
#                 logger.error(err)
#                 raise ValueError(err)
#
#     # 3. Check if workflow already in DB under a different name
#     logger.execution("Check if workflow already in DB under a different name")
#     #    Stop unless allowDuplicate flag exist
#     DBDispersionWorkflow = tk.getHermesWorkflowFromDB(dispersionWorkFlow, dispersionFlowFieldName=dispersionFlowFieldName)
#
#     if DBDispersionWorkflow is not None:
#         logger.debug(f"Found the workflow in the DB under the name {DBDispersionWorkflow[0]['desc']['workflowName']}")
#         if currentDispersionName != DBDispersionWorkflow[0]['desc']['workflowName'] and not arguments.allowDuplicate:
#             info = f"Current simulation name is {currentDispersionName}. The workflow already in DB under the name {DBDispersionWorkflow[0]['desc']['workflowName']}. use the --allowDuplicate to continue"
#             logger.error(info)
#             raise ValueError(info)
#
#     dispersionName = arguments.dispersionWorkflow.split(".")[0]
#     dispersionDirectoryName = os.path.abspath(dispersionName)
#     logger.debug(f"Getting the dispersion directory name from {arguments.dispersionWorkflow}: using {dispersionDirectoryName}")
#
#     # 4.  check if the dispersion directory exists.
#     #    if it exist, remove if the --rewrite flag exists.
#     logger.execution("Check if dispersion already exists on the disk. If it is remove for fresh start only in rewrite flag exists")
#     if os.path.isdir(dispersionDirectoryName):
#         if needRewrite:
#             logger.execution(f"The directory {dispersionDirectoryName} exists and differ than DB. Rewriting is needed as disk version was selected with --updateDB flag")
#         if arguments.rewrite:
#             logger.info(f"rewrite flag exists: removing the directory {dispersionDirectoryName}")
#             shutil.rmtree(dispersionDirectoryName)
#         else:
#             err = f"The dispersion directory {dispersionDirectoryName} already exists. use --rewrite to force removing and recreating"
#             logger.error(err)
#             raise ValueError(err)
#
#     # 3. Create the dispersion case and link to the workflow.
#     logger.info(f"Creating dispersion case {dispersionDirectoryName}  and linking to {dispersionFlowFieldDirectory}")
#     tk.stochasticLagrangian.createDispersionCaseDirectory(dispersionDirectoryName,dispersionFlowDirectory=dispersionFlowFieldDirectory)
#
#     # 4. Build the workflow and execute it.
#     logger.info("Building the workflow and executing it, forcing rebuild and reexecute")
#     arguments.force = True
#     arguments.workflow = arguments.dispersionWorkflow
#
#     handler_buildExecute(arguments)
#
#     # 5. Add/update to the DB.
#     logger.execution("add to DB if not exist. If exists, check what needs to be updated (dispersionFlowField or the entire code)")
#     if needDBAdd:
#         groupName = dispersionName.split("_")[0]
#         groupID   = dispersionName.split("_")[1]
#         logger.debug(f"Adding new document with group {groupName} and {groupID}")
#
#         parameters = DispersionWorkflow.parametersJSON
#
#         if len(DBDispersionWorkflowDocument) > 0:
#             parameters['flowParameters'] = DBDispersionWorkflowDocument[0].desc['flowParameters']
#
#         doc = tk.addSimulationsDocument(resource=dispersionFlowFieldDirectory,
#                                           dataFormat=datatypes.STRING,
#                                           type=workflowsTypes.OF_DISPERSION.value,
#                                           desc=dict(
#                                               dispersionFlowFieldName=dispersionFlowFieldName,
#                                               groupName=groupName,
#                                               groupID=groupID,
#                                               workflowName=dispersionDirectoryName,
#                                               workflowType=DispersionWorkflow.workflowType,
#                                               workflow=DispersionWorkflow.json,
#                                               parameters=parameters)
#                                           )
#
#     elif (updateDispersionFlowField or updateWorkflow):
#         doc = dispersionFlowFieldDocumentList[0]
#         if updateDispersionFlowField:
#             logger.debug("Updating dispersion flow field")
#             doc.desc['dispersionFlowFieldName'] = dispersionFlowFieldName
#
#             if len(DBDispersionWorkflowDocument) > 0:
#                 doc.desc['parameters']['flowParameters'] = DBDispersionWorkflowDocument[0].desc['flowParameters']
#
#
#         if updateWorkflow:
#             logger.debug("Updating workflow")
#             doc.desc['workflow'] = DispersionWorkflow.json
#             doc.desc['parameters'] = DispersionWorkflow.parametersJSON
#
#             if len(DBDispersionWorkflowDocument) > 0:
#                 doc.desc['parameters']['flowParameters'] = DBDispersionWorkflowDocument[0].desc['flowParameters']
#
#         logger.debug("Save to DB")
#         doc.save()
