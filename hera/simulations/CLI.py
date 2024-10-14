import json
import logging
import os

from .. import toolkitHome
from .hermesWorkflowToolkit import actionModes
from ..utils import loadJSON,compareJSONS
from hermes import workflow
from hermes.utils.workflowAssembly import handler_build,handler_buildExecute,handler_expand,handler_execute


def WorkflowsGroup_list(args):
    """
        List all the simulations group
    Parameters
    ----------
    args
        projectName : the project name. If not supplied, get from caseConfiguration

    Returns
    -------

    """
    logger = logging.getLogger("hera.bin.hera_workflows.add")
    logger.info(" -- Starting: adding workflow to the group --")

    projectName = args.projectName

    if projectName is None:
        logger.info(f"Project name is not found. Trying to load from caseConfiguration ")
        cse = loadJSON("caseConfiguration.json")
        projectName = cse.get('projectName')
        if projectName is None:
            raise ValueError(f"Project name not supplied. Provide projectName using --projectName, or in the caseConfiguration.json, projectName key ")
    else:
        logger.info(f"Using project name {projectName}. ")


    solver = args.solver
    workflowName = args.workflowName

    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=projectName)
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
        execution   : bool [default false]
            If True, then execute the workflow.

    Returns
    -------
        None
    """
    logger = logging.getLogger("hera.bin.hera_workflows.add")
    logger.info(" -- Starting: adding workflow to the group --")

    projectName = args.projectName

    if projectName is None:
        logger.info(f"Project name is not found. Trying to load from caseConfiguration ")
        cse = loadJSON("caseConfiguration.json")
        projectName = cse.get('projectName')
        if projectName is None:
            raise ValueError(f"Project name not supplied. Provide projectName using --projectName, or in the caseConfiguration.json, projectName key ")
    else:
        logger.info(f"Using project name {projectName}. ")

    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=projectName)

    workflowFile = args.workflow
    logger.info(f"Adding workflow in {workflowFile}  to DB")

    execute = args.execute
    try:
        wftk.addWorkflowToGroup(workflowJSON=workflowFile,
                                groupName= args.workflowGroup,
                                assignName=args.assignName,
                                overwrite=args.overwrite,
                                buildExecute=execute,
                                force=args.force)

    except FileExistsError as e:
        err = f"{str(e)}, use --force if you want to have duplicate records"
        logger.error(err)
        print(err)


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

    if arguments.projectName is None:
        logger.debug(f"projectName is not provided. Looking for the project name in the caseConfiguration.json file (projectName key) ")
        caseConfiguration = loadJSON("caseConfiguration.json")
        projectName = caseConfiguration['projectName']

    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=projectName)

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
    logger = logging.getLogger("hera.bin.hera_workflows.delete")
    logger.info(f" -- Starting: Deleting workflows --")


    if arguments.projectName is None:
        logger.debug(f"projectName is not provided. Looking for the project name in the caseConfiguration.json file (projectName key) ")
        caseConfiguration = loadJSON("caseConfiguration.json")
        projectName = caseConfiguration['projectName']

    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=projectName)
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
    logger = logging.getLogger("hera.bin.hera_workflows.delete")
    logger.info(f" -- Starting: Deleting workflows --")

    if arguments.projectName is None:
        logger.debug(
            f"projectName is not provided. Looking for the project name in the caseConfiguration.json file (projectName key) ")
        caseConfiguration = loadJSON("caseConfiguration.json")
        projectName = caseConfiguration['projectName']

    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=projectName)
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
            if len(res.columns)==1:
                print("\t\t ** Disk and DB are identical")
                print(" ")
            else:
                print(res)




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

    if arguments.projectName is None:
        logger.debug(f"projectName is not provided. Looking for the project name in the caseConfiguration.json file (projectName key) ")
        caseConfiguration = loadJSON("caseConfiguration.json")
        projectName = caseConfiguration['projectName']

    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=projectName)

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
        print(f"{arguments.object} is not a simulation, directory, workflow file or a simulation group in project {projectName} ")
    workflowGroup = simDocument[0].desc[wftk.DESC_GROUPNAME]

    listNodes     = arguments.nodes
    parameters    = arguments.parameters

    simulationList = wftk.listWorkflows(workflowGroup=workflowGroup,
                                        listNodes=listNodes,
                                        listParameters=parameters)


    title = f"The simulations in group *{workflowGroup}*  in project *{projectName}* "
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
        if arguments.projectName is None:
            raise ValueError("Must supply a project name for a non-file workflow")

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
    logger = logging.getLogger("hera.bin.hera_workflows.workflow_compare")
    logger.info(f" -- Starting: Listing workflow nodes --")

    if arguments.projectName is None:
        logger.debug(f"projectName is not provided. Looking for the project name in the caseConfiguration.json file (projectName key) ")
        caseConfiguration = loadJSON("caseConfiguration.json")
        projectName = caseConfiguration['projectName']

    wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName=projectName)

    res = wftk.compareWorkflow(arguments.workflows, longFormat=arguments.longFormat, transpose=arguments.transpose)


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
        print(f"Could not found any workflows to compare in project {projectName}")
    else:
        print(output)

        if arguments.file is not None:
            flName = arguments.file if "." in arguments.file else f"{arguments.file}.{ext}"

            with open(flName,"w") as outputFile:
                outputFile.write(output)

def workflow_expand(arguments):
    parser_expandWorkflow(arguments)

def workflow_build(arguments):
    handler_build(arguments)

def workflow_execute(arguments):
    handler_execute(arguments)

def workflow_buildExecute(arguments):
    handler_buildExecute(arguments)

