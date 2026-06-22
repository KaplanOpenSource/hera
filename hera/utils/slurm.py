from hera.utils.logging import get_logger

def prepareSlurmScriptExecution(*, scriptPath=None, script=None,
                            slurmExecutionFilePath="submit_all.sh",
                            jobDirListFilePath="jobPaths.txt",
                            allocateProcessorsPerRun=None,
                            memoryInGB=None,
                            jobName="generic_job",
                            quiet=False,
                            exclusive=False):
    """
        Prepare slurm batch file to run a script provided either as a file or str on all jobs in jobListFile.
        Note that script can access all parameters of the SBATCH file including $dir an {{SLURM_ARRAY_TASK_ID}}

        This assumes each job is ran from a specific folder for that job and cd's to that folder prior

    Parameters
    ----------
    scriptPath : str | None
            Path to a script to run per job in jobDirListFilePath
    script : str | None
            Script to run per job in jobDirListFilePath
    slurmExecutionFileName: str
            The name of the bash file to create with the slurm batch run
    jobDirListFilePath:
            path to file with all the paths to the directories for each job separated with new lines
    allocateProcessorsPerRun : int | None
            How many nodes(currently just threads) are used per job, in slurm documentation they describe "This option advises ... that job steps run ... will launch a maximum of number tasks and to provide for sufficient resources."
    memoryInGB : str | int | None
            Should memory be limited, affect depends on configuration(might kill if memory use is exceeded)
    jobName : bool
            name for slurm job
    exclusive : bool
            Should slurm run one job at a time on a GRES(Generic RESource), means server/computer

    Returns
    -------

    """

    logger = get_logger(None,"prepareSlurmScriptExecution")
    if scriptPath is None and script is None:
        logger.error("either scriptPath or script must be provided")
        return
    elif not(scriptPath is None) and not(script is None):
        logger.error("only one of scriptPath or script must be provided")
        return
    if scriptPath is not None:
        with open(scriptPath, 'r') as scriptFile:
            script = scriptFile.read()
            
    try:
        with open(jobDirListFilePath, 'r') as file:
            counter = sum(1 for _ in file)
    except FileNotFoundError:
        logger.error(f"case list file({jobDirListFilePath}) not found")
        return
    except IOError:
        logger.error(f"unable to read case list file({jobDirListFilePath})")
        return

    processorLine = f"#SBATCH -n {allocateProcessorsPerRun} # number of CPUs per job" if allocateProcessorsPerRun is not None else ""
    meomoryLine   = f"#SBATCH -mem={memoryInGB}G" if memoryInGB is not None else ""
    exclusiveLine = f"#SBATCH --exclusive" if exclusive is not None and exclusive else ""
    slurmBatchFile = f"""#!/bin/bash
#SBATCH --job-name={jobName}
#SBATCH --array=1-{counter}     
{processorLine}
{meomoryLine}
{exclusiveLine}
#SBATCH --output=slurm_%A_%a.out

# Read the directory for this array task
dir=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {jobDirListFilePath})

{'echo "Running script for job in directory: $dir"' if not quiet else ""}

{script}
"""
    logger.info(f"Writing slurm execution submit file to {slurmExecutionFilePath}")
    with open(slurmExecutionFilePath,"w") as outputFile:
        outputFile.write(slurmBatchFile)