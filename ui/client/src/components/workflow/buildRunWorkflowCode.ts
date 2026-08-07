// Python that builds and executes a saved Hermes workflow from the DB, the same
// as the CLI `hera-workflows buildExecute`. Runs with the local Luigi scheduler
// (no separate process). Returns the dispatch id of the run.
export const buildRunWorkflowCode = ({
  projectName,
  workflowName,
}: {
  projectName: string,
  workflowName: string,
}) => {
  // Temporarily add the toolkit's files directory to PYTHONPATH (restored after) so the Luigi subprocess can import the generated workflow module.
  return `
import os
from hera import toolkitHome
wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName='${projectName}')
_old_pythonpath = os.environ.get('PYTHONPATH')
os.environ['PYTHONPATH'] = wftk.FilesDirectory + os.pathsep + (_old_pythonpath or '')
try:
    dispatch_id = wftk.executeWorkflowFromDB('${workflowName}', scheduler='local')
finally:
    if _old_pythonpath is None:
        os.environ.pop('PYTHONPATH', None)
    else:
        os.environ['PYTHONPATH'] = _old_pythonpath
`;
};
