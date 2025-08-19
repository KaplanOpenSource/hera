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
from .preprocessOFObjects import OFObjectHome
from ..CLI  import workflow_add


def Foam_createEmpty(arguments):
    logger = logging.getLogger("hera.bin")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    if 'projectName' not in arguments:
        configurationFile = arguments.configurationFile if 'configurationFile' in arguments else "caseConfiguration.json"
        try:
            configuration = loadJSON(configurationFile)
        except:
            err = f"Configuration file {configurationFile} not found! creating a basic file"
            with open(configurationFile, 'w') as defaultConfFile:
                json.dump(dict(projectName=None), defaultConfFile, indent=4)

            raise ValueError(err)

        projectName = configuration['projectName']
    else:
        projectName = arguments.projectName

    logger.info(f"Using project {projectName}")
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)

    simulationType = tk.FLOWTYPE_INCOMPRESSIBLE if arguments.incompressible else tk.FLOWTYPE_COMPRESSIBLE

    tk.createEmptyCase(caseDirectory=arguments.caseDirectory,
                       fieldList=arguments.fields,
                       flowType=simulationType,
                       additionalFieldsDescription=arguments.fieldsDescription)

    logger.debug(f"----- End -----")


def Foam_parser_FieldDescription(arguments):
    logger = logging.getLogger("hera.bin")
    logger.debug(f"----- Start : Foam_parser_FieldDescription -----")
    logger.debug(f" arguments: {arguments}")

    if arguments.fields is not None:
        field_len = len(arguments.fields)
    else:
        field_len = 0


    jsonExample = dict()
    if field_len > 0:
        for fieldName in arguments.fields:
            jsonExample[fieldName] = dict(dimensions=OFObjectHome.getDimensions(), componentNames=None)
    else:
        jsonExample["exampleField"] = dict(dimensions=OFObjectHome.getDimensions(), componentNames=None)

    with open(arguments.fileName, "w") as outfile:
        json.dump(jsonExample, outfile, indent=4)

def foam_solver_template_buildExecute(arguments):
    """
        Adds the workflow to the DB and executes it.

        Can be invoked without executing or without adding to the db.

    Parameters
    ----------
    arguments

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin")
    logger.debug(f"----- Start : foam_templates_flow_list-----")
    logger.debug(f" arguments: {arguments}")

    arguments.force = True
    if arguments.noDB:
        from hermes.utils.workflowAssembly import handler_buildExecute
        handler_buildExecute(arguments)
    else:
        workflow_add(arguments)
        handler_buildExecute(arguments)

def foam_solver_templates_list(arguments):
    logger = logging.getLogger("hera.bin")
    logger.debug(f"----- Start : foam_templates_flow_list-----")
    logger.debug(f" arguments: {arguments}")

    projectName = None if 'projectName' not in arguments else arguments.projectName  # from hera 2.13.2 the toolkit searches the project name in the case file.
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    templates = tk.listHermesSolverTemplates(arguments.solver)

    ttl = f"The templates for project {tk.projectName} with solver {arguments.solver}"
    print()
    print("-" * len(ttl))
    print(ttl)
    print("-" * len(ttl))
    print(templates)


def foam_solver_template_create(arguments):
    logger = logging.getLogger("hera.bin")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    if 'projectName' in arguments:
        logger.debug("Take the projectName from the arguments")
        projectName = arguments.projectName
    else:
        logger.debug("Take the projectName from the directory")
        projectName = None

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    templates = tk.listHermesSolverTemplates(arguments.solver)
    if arguments.templateName not in templates.index:
        err = f"{arguments.templateName} is not known. Use one of the \n " + str(templates)
        logger.error(err)
        raise ValueError(err)

    outputPath = os.getcwd() if arguments.projectPath is None else arguments.projectPath
    groupName = arguments.templateName if arguments.groupName is None else arguments.groupName
    logger.info(f"Saving the simulation with the group name: {groupName}")

    tname = arguments.templateName
    dataSourceName = templates.query(f"name==@tname").iloc[0].datasourceName
    with open(os.path.join(outputPath, f"{groupName}_1.json"), "w") as outFile:
        json.dump(tk.getDataSourceData(dataSourceName), outFile, indent=4)


def foam_templates_node_list(arguments):
    logger = logging.getLogger("hera.bin")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    projectName = None if 'projectName' not in arguments else arguments.projectName  # from hera 2.13.2 the toolkit searches the project name in the case file.
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    templates = tk.listHermesNodesTemplates()

    ttl = f"The templates of nodes for project {tk.projectName} with solver {arguments.solver}"
    print()
    print("-" * len(ttl))
    print(ttl)
    print("-" * len(ttl))
    print(templates)

def foam_solver_simulations_list(arguments):
    logger = logging.getLogger("hera.bin")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    projectName = None if 'projectName' not in arguments else arguments.projectName  # from hera 2.13.2 the toolkit searches the project name in the case file.
    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)

    simDocument = wftk.getWorkflowListOfSolvers(arguments.solver)

    if len(simDocument) == 0:
        print(f"There are no cases for {arguments.solver} in {projectName} ")

    for groupName in set([x.desc['groupName'] for x in simDocument]):
        ttl = " "*10 + f"Group name: {groupName}" + " "*10
        print(ttl)
        print(f"-"*len(ttl))

        res = wftk.compareWorkflows([groupName], longFormat=arguments.longFormat, transpose=arguments.transpose)

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


def stochasticLagrangian_dispersionFlow_create(arguments):
    """

    Parameters
    ----------
    arguments

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    projectName = arguments.projectName
    logger.info(f"Adding dispersion flow to project {projectName}. Overwriting? {arguments.overwrite}")
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=arguments.projectName)

    params = arguments.dispersionFlowParams.replace("'",'"').replace("True","true").replace("False","false").replace("None","null")
    flowdata = loadJSON(params)
    flowName = arguments.dispersionFlowName
    logger.info(f"Createing dispersion flow (DFF) {flowName} for the flow {arguments.OriginalFlowField}")

    # dispersionFieldList = []
    # if 'dispersionFlowFields' in arguments:
    #     if arguments.dispersionFlowFields is not None:
    #         dispersionFieldList = list(numpy.atleast_1d(arguments.dispersionFlowFields))

    try:
        tk.stochasticLagrangian.createDispersionFlowField(flowName=flowName,
                                                          flowData=flowdata,
                                                          OriginalFlowField=arguments.OriginalFlowField,
                                                          dispersionDuration=arguments.dispersionDuration,
                                                          overwrite=arguments.overwrite)
    except FileExistsError:
        err = f"Flow field {flowName} Already exists. Use --overwrite to recreate"
        logger.error(err)

    logger.debug(f"----- End -----")

def foam_mesh_blockMesh(arguments):
    raise NotImplementedError("Not implemented yet, blockMesh")

def foam_mesh_setDomainHeight(arguments):
    raise NotImplementedError("Not implemented yet, maybe we will use classy blocks")

def IC_hydrostaticPressure(argumets):
    raise NotImplementedError("Use the openfoam application")

def stochasticLagrangian_dispersionFlow_writeEmptyTemplate(arguments):
    """
        Writes an empty JSON file for the disprsion flow.
        The name of the dispersion workflow is the file name

        Also allows the setting of a constant value to the variable.
        For more compllicated assigments use setFields, or the python interface of hera.

    Parameters
    ----------
    arguments

    Returns
    -------

    """
    flowData = {
        "originalFlow": {
            "time": {
                "temporalType": "steadyState|dynamic",
                "timestep": "< time >"
            },
            "linkMeshSymbolically": True
        },
        "dispersionDuration": "< duration >",
        "dispersionFields" : {
            "<field name>" : "<constant value>"
        }
    }
    with open(arguments.templateFile,'w') as outFile:
        json.dump(flowData,outFile,indent=4)

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
    logger.debug(f"----- Start -----")
    logger.debug(f" makeDispersionFlow with arguments: {arguments}")

    if 'projectName' not in arguments:
        configurationFile = arguments.configurationFile if 'configurationFile' in arguments else "caseConfiguration.json"

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
    logger.debug(f"----- Start -----")
    logger.debug(f"  Got arguments: {arguments}")

    if ('projectName' in arguments) and (arguments.projectName is not None):
        projectName = arguments.projectName
    else:
        configurationFile = arguments.configurationFile if 'configurationFile' in arguments else "caseConfiguration.json"
        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']

    logger.info(f"Using project {projectName}")

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)

    dispersionDirectoryName = arguments.dispersionName
    logger.info(f"Creating dispersion case {dispersionDirectoryName}")

    dispersionFlowFieldName = arguments.dispersionFlowField
    logger.info(f"Getting the dispersion flowField {dispersionFlowFieldName} ")
    doc = tk.getWorkflowDocumentFromDB(dispersionFlowFieldName, tk.DOCTYPE_OF_FLOWDISPERSION)
    if len(doc) == 0 or arguments.overwrite:
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
            shutil.rmtree(dispersionDirectoryName)
        else:
            err = "Dispersion flow exists. Use --overwrite to force recreation"
            raise ValueError(err)

    # 3. Create the dispersion case and link to the workflow.
    logger.info(f"Creating dispersion case {dispersionDirectoryName}  and linking to {dispersionFlowField}")
    tk.stochasticLagrangian.createAndLinkDispersionCaseDirectory(dispersionDirectoryName,
                                                                 dispersionFlowDirectory=dispersionFlowField)

def stochasticLagrangian_source_cylinder(arguments):
    logger = logging.getLogger("hera.bin.stochasticLagrangian_source_cylinder")
    logger.debug(f"----- Start -----")
    logger.debug(f"  Got arguments: {arguments}")

    center = arguments.center
    params = dict(x=float(center[0]),
                  y=float(center[1]),
                  z=float(center[2]),
                  radius=float(arguments.radius),
                  height=float(arguments.height),
                  nParticles=int(arguments.particles)
                  )
    if ('projectName' in arguments) and (arguments.projectName is not None):
        projectName = arguments.projectName
    else:
        configurationFile = arguments.configurationFile if 'configurationFile' in arguments else "caseConfiguration.json"
        logger.debug(f"Loading configuration file {configurationFile}")
        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']

    dispersionName = arguments.dispersionName

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    tk.stochasticLagrangian.writeParticlePositionFile(type="Cylinder", dispersionName=dispersionName, **params)

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

    case = os.path.abspath(args.casePath)
    massFileName = f"{args.patch}Mass" if args.massFileName is None else args.massFileName
    dt = args.dt
    LSMtoolkit = toolkitHome.getToolkit(toolkitHome.OF_LSM, "tmpProject", casePath=case)
    data = LSMtoolkit.analysis.getMassFromLog(logFile=args.logFile, solver=args.solver)
    data = data.loc[data.filterType == args.patch].loc[data.action == args.action]
    data["diffMass"] = data.mass.diff()
    data = data.fillna(0)
    if dt is None:
        timesteps = data.time
    else:
        timesteps = numpy.arange(data.time.min(), data.time.max(), float(dt))
        times = pandas.DataFrame({"time": timesteps})
        data = data.set_index("time").join(times.set_index("time"), how="outer").reset_index()
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
        newstr += f"{float(data.loc[data.time == time].diffMass)}\n"
    newstr += ")"
    print("saving in ", os.path.join(case, "constant", massFileName))
    f = open(os.path.join(case, "constant", massFileName), "w")
    f.write(newstr)
    f.close()

def stochasticLagrangian_postProcess_toParquet(arguments):
    """
        Creates a parquet file from the case results
    """
    logger = logging.getLogger("hera.bin.stochasticLagrangian_postProcess_toParquet")
    logger.debug(f"----- Start -----")
    logger.debug(f"  Got arguments: {arguments}")

    if ('projectName' in arguments) and (arguments.projectName is not None):
        projectName = arguments.projectName
    else:
        configurationFile = arguments.configurationFile if 'configurationFile' in arguments else "caseConfiguration.json"
        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']

    logger.info(f"Using project {projectName} for cloud name {arguments.cloudName}")

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    docList = tk.getWorkflowDocumentFromDB(arguments.dispersionName)
    if len(docList) == 0:
        logger.info(f"The name {arguments.dispersionName} is not found. Try to use it as a directory")
        if os.path.isdir(arguments.dispersionName):
            logger.info(
                f"Found {arguments.dispersionName} as directory. Trying to use it as lagrangian simulation and save it in VTK format")
            outputName = arguments.dispersionName
            cache = False
    else:
        doc = docList[0]
        outputName = doc.desc['workflowName']
        cache = True



    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    data = tk.stochasticLagrangian.getCaseResults(caseDescriptor=outputName, withVelocity=True, withMass=True,
                                                  overwrite=arguments.overwrite, cloudName=arguments.cloudName)


def stochasticLagrangian_postProcess_toVTK(arguments):
    """
        Creates VTK files from the case.
        Creates the parquet file if does not exist.
    """
    logger = logging.getLogger("hera.bin.stochasticLagrangian_postProcess_toParquet")
    logger.debug(f"----- Start -----")
    logger.debug(f"  Got arguments: {arguments}")

    if ('projectName' in arguments) and (arguments.projectName is not None):
        projectName = arguments.projectName
    else:
        configurationFile = arguments.configurationFile if 'configurationFile' in arguments else "caseConfiguration.json"
        configuration = loadJSON(configurationFile)
        projectName = configuration['projectName']

    logger.info(f"Using project {projectName}")

    if 'outputdir' in arguments:
        base_outputdir = os.getcwd() if arguments.outputdir is None else arguments.outputdir
    else:
        base_outputdir = os.getcwd()

    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
    docList = tk.getWorkflowDocumentFromDB(arguments.dispersionName)
    if len(docList) == 0:
        logger.info(f"The name {arguments.dispersionName} is not found. Try to use it as a directory")
        if os.path.isdir(arguments.dispersionName):
            logger.info(
                f"Found {arguments.dispersionName} as directory. Trying to use it as lagrangian simulation and save it in VTK format")
            outputName = arguments.dispersionName
    else:
        doc = docList[0]
        outputName = doc.desc['workflowName']
        logger.info(f"The name {arguments.dispersionName} found. Use the workflow name {outputName}")

    outputdir = os.path.join(base_outputdir, "VTK", outputName)
    logger.info(f"Writing values to directory {outputdir}")
    os.makedirs(outputdir, exist_ok=True)
    data = tk.stochasticLagrangian.getCaseResults(caseDescriptor=outputName, withVelocity=True, withMass=True,overwrite=arguments.overwrite,cache=False)
    tk.presentation.toUnstructuredVTK(data=data, outputdirectory=outputdir, filename=arguments.cloudName,
                                      overwrite=arguments.overwrite, xcoord="x", ycoord="y", zcoord="z")

def objects_createVerticesAndBoundary(arguments):
    logger = logging.getLogger("hera.bin")
    logger.debug(f"----- Start -----")
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
    objFile = FreeCAD.getDocument("Unnamed")

    ret = dict()
    for fieldName in arguments.fields:
        boundary = dict()
        for regionObj in objFile.findObjects():
            boundary[regionObj.Name] = dict(type="zeroGradient")
        ret[fieldName] = dict(boundaryField=boundary)

    print("----- Bounding box vertices -----------")
    corners = getObjFileBoundaries(fileName)

    xList = ['XMin', 'XMax', 'XMax', 'XMin', 'XMin', 'XMax', 'XMax', 'XMin']
    yList = ['YMin', 'YMin', 'YMax', 'YMax', 'YMin', 'YMin', 'YMax', 'YMax']
    zList = ['ZMin', 'ZMin', 'ZMin', 'ZMin', 'ZMax', 'ZMax', 'ZMax', 'ZMax']

    verticesList = []
    for x, y, z in zip(xList, yList, zList):
        xshift = 0.1 if 'Max' in x else -0.1
        yshift = 0.1 if 'Max' in y else -0.1
        zshift = 0.1 if 'Max' in z else -0.1

        verticesList.append([corners[x] + xshift, corners[y] + yshift, corners[z] + zshift])
    print(json.dumps(dict(vertices=verticesList), indent=4))

    print("\n\n\n")
    print("----- Boundary conditions -------------")
    print(json.dumps(ret, indent=4, sort_keys=True))

############################################################## Workflow.
def foam_mesh_blockMesh(arguments):
    """
        Adjusts the blockMesh according to the stl file.

    Parameters
    ----------
    arguments

    Returns
    -------

    """
    pass

def foam_mesh_setDomainHeight(arguments):
    """
        Sets the height of the domain and updates the number of cells.
    Parameters
    ----------
    arguments

    Returns
    -------

    """
    pass

def foam_snappyhexmesh_addobject(arguments):
    """
        Adding the snappy hex mesh node (if does not exist).
        If exists, just adds the object.

    Parameters
    ----------
    arguments

    Returns
    -------

    """
    pass

def foam_snappyhexmesh_setLocationInDomain(arguments):
    """
        Sets the location in mesh.
    Parameters
    ----------
    arguments

    Returns
    -------

    """
    pass

def foam_IC(arguments):
    """
        Sets the Initial conditions.
    Parameters
    ----------
    arguments

    Returns
    -------

    """
    pass

def IC_hydrostaticPressure(arguments):
    """
        Sets the hydrostatic pressure for the p variable.
    Parameters
    ----------
    arguments

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin")
    logger.debug(f"----- Start -----")
    logger.debug(f" arguments: {arguments}")

    if 'projectName' not in arguments:
        configurationFile = arguments.configurationFile if 'configurationFile' in arguments else "caseConfiguration.json"
        try:
            configuration = loadJSON(configurationFile)
        except:
            err = f"Configuration file {configurationFile} not found! creating a basic file"
            with open(configurationFile, 'w') as defaultConfFile:
                json.dump(dict(projectName=None), defaultConfFile, indent=4)

            raise ValueError(err)

        projectName = configuration['projectName']
    else:
        projectName = arguments.projectName

    logger.info(f"Using project {projectName}")
    tk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)

    simulationType = tk.FLOWTYPE_INCOMPRESSIBLE if arguments.incompressible else tk.FLOWTYPE_COMPRESSIBLE

    cellCenters = tk.getMesh(arguments.caseDirectory)
    pField = getattr(tk, arguments.solver).IC_getHydrostaticPressure(arguments.caseDirectory)
    pField.writeToCase(arguments.caseDirectory, timeOrLocation=arguments.startTime)

def foam_BC(arguments):
    """
        Sets the boundary codntions.
    Parameters
    ----------
    arguments

    Returns
    -------

    """
    pass

def hermes_buildExecute(arguments):
    logger = logging.getLogger("bin.hermes_buildExecute")
    try:
        from hermes.utils.workflowAssembly import handler_build, handler_buildExecute, handler_expand, \
            handler_execute

    except ImportError as e:
        err = f"hermes package is not found. Cannot build and execute file. "
        logger.error(err)
        raise ImportError(err)

    handler_buildExecute(arguments)

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
#     logger.debug(f"----- Start -----")
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
#         logger.debug(f" the workflow {currentDispersionName} in the DB. First see if the dispersion flow field is similar")
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
#         logger.debug(f"Check if the disk workflow is identical to the DB workflow")
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
#                 logger.debug(f"Exporting the workflow to file {arguments.dispersionWorkflow}")
#                 with open(arguments.dispersionWorkflow, 'w') as JSONOut:
#                     json.dump(DBDispersionWorkflowDocument[0]['desc']['workflow'],JSONOut,indent=4)
#
#                 DispersionWorkflow = DBDispersionWorkflow
#             elif arguments.updateDB:
#                 logger.debug("Updating the DB with the workflow from disk")
#                 DispersionWorkflow = diskDispersionWorkflow
#                 updateWorkflow = True
#                 needRewrite = True
#             else:
#                 err = f"The workflow in {arguments.dispersionWorkflow} file and in the DB are different. Use the --updateDB to update the DB or --exportFromDB to rewrite the file"
#                 logger.error(err)
#                 raise ValueError(err)
#
#     # 3. Check if workflow already in DB under a different name
#     logger.debug("Check if workflow already in DB under a different name")
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
#     logger.debug("Check if dispersion already exists on the disk. If it is remove for fresh start only in rewrite flag exists")
#     if os.path.isdir(dispersionDirectoryName):
#         if needRewrite:
#             logger.debug(f"The directory {dispersionDirectoryName} exists and differ than DB. Rewriting is needed as disk version was selected with --updateDB flag")
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
#     logger.debug("add to DB if not exist. If exists, check what needs to be updated (dispersionFlowField or the entire code)")
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
