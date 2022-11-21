# import pandas
import numpy
import os
from . import HERAMETADATA
from ...utils import loggedObject,loadJSON
# from ..utils.coordinateHandler import coordinateHandler
# handler = coordinateHandler()


def getNumberOfSubdomains(caseDirectory):
    """
        Reads the decomposeParDict and returns the number of subdomains of that particular case.
    Parameters
    ----------
    caseDirectory

    Returns
    -------

    """
    decomposeParDictFileName = os.join(caseDirectory,"system","decomposeParDict")
    from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

    f=ParsedParameterFile(decomposeParDictFileName)
    return f['numberOfSubdomains']

def buildCaseExecutionScript(caseDirectory, workflow, isParallel=None, isSlurm=None):
    """
        Writes the allRun file that executes the workflow and the allClean file
        that cleans the case to the case directory.

    Parameters
    ----------
    caseDirectory: str
        The directory to write the files to.

    workflow: dict
        A configuration file that is used to build the execution of the node.
        Specifically, includes a node 'caseExecution' with the structure:

         {
            "parallelCase": true|false             # Write the execution for parallel execution.
            "runFile": [
                        --------------- A list of nodes.
              {
                "name": "blockMesh",               # The name of the program to execute.
                "couldRunInParallel": false,       # Write as parallel (only if parallel case is True).
                "parameters": null                 # Parameters for each run.
              }
            ]
                        --------------- A list of nodes.
         }

    isParallel  : bool
            If true, force parallel, if False force single.
            If None, get from file.

    isSlurm : bool
            If true, use SLURM suuport to allocate and run simulation.


    """
    logger = loggedObject(loggerName="simulations.openFoam.createRunScripts").logger

    logger.info("---- Start ----")


    workflow = loadJSON(workflow)
    try:
        execConfiguration = workflow[HERAMETADATA]['caseExecution']
    except KeyError:
        err = "'caseExecution' node not found in configuration"
        logger.error(err)
        raise ValueError(err)

    isSlurm   = execConfiguration.get('slurm',False) if isSlurm is None else isSlurm

    isParallel = execConfiguration['parallelCase'] if isParallel is None else isParallel
    logger.info(f"Running this case as parallel? {isParallel}")

    execLine = ""

    for execNode in execConfiguration['runFile']:
        logger.execution(f"Processing Node {execNode['name']}")

        parallelFlag = "-parallel" if (isParallel and execNode['couldRunInParallel']) else ""
        progName = execNode['name']
        parameters = execNode.get('parameters',None)

        if parameters is not None:
            params   = " ".join(numpy.atleast_1d(execNode['parameters']))
        else:
            params = ""

        foamJob = execNode.get("foamJob",True)
        slurm   = ""
        if isSlurm:
            slurm = "srun"
            if isParallel:
                procCount = getNumberOfSubdomains(caseDirectory)
                execLine += f"salloc {procCount}\n"

        if foamJob:
            execLine += f"{slurm} foamJob {parallelFlag} -append -screen -wait {progName} {params}\n"
        else:
            execLine += f"{slurm} {progName} {params}\n"

    allrunFile = os.path.join(caseDirectory,"Allrun")
    with open(allrunFile,'w') as execFile:
        execFile.write(execLine)
    os.chmod(allrunFile, 0o777)

    # Now write the allClean file.
    allCleanContent = """
#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

cp 0.parallel/* 0
cleanCase
    """
    allcleanFile = os.path.join(caseDirectory,"Allclean")
    with open(allcleanFile,'w') as allclean:
        allclean.write(allCleanContent)

    os.chmod(allcleanFile, 0o777)


def getCellDataAndGroundData(casePath,ground="ground"):
    f = open(os.path.join(casePath, "0", "cellCenters"), "r")
    lines = f.readlines()
    f.close()
    fboundary = open(os.path.join(casePath, "0", "Hmix"), "r")
    boundarylines = fboundary.readlines()
    fboundary.close()

    for i in range(len(lines)):
        if "internalField" in lines[i]:
            cellStart = i + 3
            nCells = int(lines[i + 1])
        if ground in lines[i]:
            break
    for boundaryLine in range(len(boundarylines)):
        if "boundaryField" in boundarylines[boundaryLine]:
            break

    nGroundValues = int(lines[i + 4])
    groundData = pandas.read_csv(os.path.join(casePath, "0", "cellCenters"), skiprows=i + 6,
                                 skipfooter=len(lines) - (i + 6 + nGroundValues),
                                 engine='python',
                                 header=None,
                                 delim_whitespace=True, names=['x', 'y', 'z'])

    groundData['x'] = groundData['x'].str[1:]
    groundData['z'] = groundData['z'].str[:-1]
    groundData = groundData.astype(float)

    cellData = pandas.read_csv(os.path.join(casePath, "0", "cellCenters"), skiprows=cellStart,
                               skipfooter=len(lines) - (cellStart + nCells),
                               engine='python',
                               header=None,
                               delim_whitespace=True, names=['x', 'y', 'z'])
    cellData['x'] = cellData['x'].str[1:]
    cellData['z'] = cellData['z'].str[:-1]
    cellData = cellData.astype(float)
    return cellData,groundData