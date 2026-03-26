# Slurm Batch Jobs

The `hera.utils.slurm` module provides `prepareSlurmScriptExecution`, a function that generates a Slurm array-job submission script. Each array task runs a user-provided script from a different working directory, making it straightforward to launch many independent simulations.

## Importing

```python
from hera.utils.slurm import prepareSlurmScriptExecution
```

## Basic Usage

### 1. Prepare a Job Directory List

Create a text file with one directory path per line. Each line corresponds to one array task:

```text
/data/runs/case_001
/data/runs/case_002
/data/runs/case_003
```

### 2. Generate the Submission Script

```python
from hera.utils.slurm import prepareSlurmScriptExecution

prepareSlurmScriptExecution(
    script="cd $dir && python run_simulation.py",
    jobDirListFilePath="jobPaths.txt",
    slurmExecutionFilePath="submit_all.sh",
    jobName="my_simulations",
)
```

This writes `submit_all.sh` with Slurm directives for a 3-task array job.

### 3. Submit

```bash
sbatch submit_all.sh
```

## Using a Script File

Instead of an inline script string, point to a shell script:

```python
prepareSlurmScriptExecution(
    scriptPath="run_case.sh",
    jobDirListFilePath="jobPaths.txt",
    slurmExecutionFilePath="submit_all.sh",
    jobName="batch_run",
)
```

The contents of `run_case.sh` are embedded directly into the generated Slurm script. The variable `$dir` is available inside the script and points to the directory for the current array task.

## Parameters

| Parameter | Type | Description |
|---|---|---|
| `script` | `str` | Inline script to run per job. Mutually exclusive with `scriptPath`. |
| `scriptPath` | `str` | Path to a script file. Mutually exclusive with `script`. |
| `slurmExecutionFilePath` | `str` | Output path for the generated submission script (default: `"submit_all.sh"`). |
| `jobDirListFilePath` | `str` | Path to a text file listing one job directory per line (default: `"jobPaths.txt"`). |
| `jobName` | `str` | Name shown in `squeue` output (default: `"generic_job"`). |
| `allocateProcessorsPerRun` | `int` | Number of CPUs per array task. Omit to use the cluster default. |
| `memoryInGB` | `int` | Memory limit in GB per task. Omit for no limit. |
| `exclusive` | `bool` | If `True`, each task gets exclusive access to its node (default: `False`). |
| `quiet` | `bool` | If `False`, the script echoes the current directory before running (default: `False`). |

## Full Example with Resource Constraints

```python
prepareSlurmScriptExecution(
    script="""
cd $dir
source activate myenv
python run.py --config config.json > output.log 2>&1
""",
    jobDirListFilePath="all_cases.txt",
    slurmExecutionFilePath="submit.sh",
    jobName="cfd_sweep",
    allocateProcessorsPerRun=8,
    memoryInGB=16,
    exclusive=True,
)
```

The generated script will include:

```bash
#SBATCH -n 8 # number of CPUs per job
#SBATCH -mem=16G
#SBATCH --exclusive
```

## Integration with LSM Toolkit

When using the LSM (Land Surface Model) toolkit, the toolkit can generate the job directory list. You then pass that file to `prepareSlurmScriptExecution` to create the submission script:

```python
# After the toolkit creates case directories and writes jobPaths.txt
prepareSlurmScriptExecution(
    scriptPath="run_lsm.sh",
    jobDirListFilePath="jobPaths.txt",
    slurmExecutionFilePath="submit_lsm.sh",
    jobName="lsm_batch",
    allocateProcessorsPerRun=4,
)
```
