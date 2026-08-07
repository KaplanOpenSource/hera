// Python that builds and executes a saved Hermes workflow from the DB, the same
// as the CLI `hera-workflows buildExecute`. Runs with the local Luigi scheduler
// (no separate process). Returns the dispatch id and the captured console output.
export const buildRunWorkflowCode = ({
  projectName,
  workflowName,
}: {
  projectName: string,
  workflowName: string,
}) => {
  // Add the files dir to PYTHONPATH (so the Luigi subprocess can import the generated module) and
  // capture stdout+stderr at the fd level (the subprocess and node os.system write to fds, not sys.stdout);
  // both are restored afterwards.
  return `
import os
import sys
import tempfile
from hera import toolkitHome

workflowToolkit = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName='${projectName}')

previous_pythonpath = os.environ.get('PYTHONPATH')
os.environ['PYTHONPATH'] = workflowToolkit.FilesDirectory + os.pathsep + (previous_pythonpath or '')

log_fd, log_path = tempfile.mkstemp(suffix='.log')
saved_stdout_fd = os.dup(1)
saved_stderr_fd = os.dup(2)
sys.stdout.flush()
sys.stderr.flush()
os.dup2(log_fd, 1)
os.dup2(log_fd, 2)

try:
    dispatch_id = workflowToolkit.executeWorkflowFromDB('${workflowName}', scheduler='local')
finally:
    sys.stdout.flush()
    sys.stderr.flush()
    os.dup2(saved_stdout_fd, 1)
    os.dup2(saved_stderr_fd, 2)
    os.close(saved_stdout_fd)
    os.close(saved_stderr_fd)
    os.close(log_fd)
    if previous_pythonpath is None:
        os.environ.pop('PYTHONPATH', None)
    else:
        os.environ['PYTHONPATH'] = previous_pythonpath

with open(log_path) as log_file:
    output = log_file.read()
os.remove(log_path)
`;
};
