# Slurm Batch Jobs

The `hera.utils.slurm` module provides `prepareSlurmScriptExecution`, a function that generates a Slurm array-job submission script. Each array task reads one directory from a list file into the shell variable `$dir` and then runs a user-provided script, making it straightforward to launch many independent simulations. The generated script does **not** change into `$dir` itself; start your script with `cd "$dir"`.

## Importing

```python
from hera.utils.slurm import prepareSlurmScriptExecution
```

## Basic Usage

### 1. Prepare a Job Directory List

Create a text file with one directory path per line. Each line corresponds to one array task (task `N` reads line `N` with `sed`), so avoid blank lines:

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

This writes `submit_all.sh` with Slurm directives for a 3-task array job (`#SBATCH --array=1-3`), with per-task output in `slurm_%A_%a.out`. The list file path is written into the script as given, so use an absolute path or run `sbatch` from the same directory.

All arguments are keyword-only. On invalid input (both or neither of `script` / `scriptPath`, or a missing list file) the function logs an error and returns `None` without writing anything.

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
| `quiet` | `bool` | If `False`, the script echoes the job directory (`$dir`) before running (default: `False`). |

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

!!! warning "Known issue: memory directive"
    The memory line is written as `#SBATCH -mem=...` (single dash). The sbatch option is `--mem`, so the memory limit is not applied as intended. Until this is fixed, omit `memoryInGB` and add `#SBATCH --mem=16G` at the top of your `script` text, or edit the generated file.

## Integration with Toolkits

Two toolkits call `prepareSlurmScriptExecution` for you, after expanding a base configuration with [`JSONVariations`](json.md):

- **LSM**: `LSMToolkit.prepareSlurmLSMExecution(baseParameters, jsonVariations, templateName, ..., slurmExecutionFileName="submit_all.sh", caseListFileName="cases.txt", allocateProcessorsPerRun=None, memoryInGB=None, jobName="lsm_cases", exclusive=False)` writes one `LSM_Simulation_<i>/run_sim.py` per variation plus `cases.txt` and `submit_all.sh` under `<filesDirectory>/simulationsScripts/`.
- **OpenFOAM**: `OFToolkit.prepareSlurmWorkflowExecution(baseConfiguration, jsonVariations, slurmExecutionFileName="submit_all.sh", caseListFileName="cases.txt", allocateProcessorsPerRun=None, memoryInGB=None, jobName="foam_cases", exclusive=False, addAllRun=True)` adds each variation as a workflow to the group and writes `cases.txt` and `submit_all.sh` in the toolkit's `filesDirectory`. Each task runs `hera-workflows sync` / `buildExecute` and, when `addAllRun=True`, `bash Allrun`.

```python
lsm = toolkitHome.getToolkit(toolkitHome.LSM, projectName="MY_PROJECT")
lsm.prepareSlurmLSMExecution(
    baseParameters=base,
    jsonVariations=variations,
    templateName="my_template",
    allocateProcessorsPerRun=4,
)
# then: cd <filesDirectory>/simulationsScripts && sbatch submit_all.sh
```
