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
  return `
from hera import toolkitHome
wftk = toolkitHome.getToolkit(toolkitName=toolkitHome.SIMULATIONS_WORKFLOWS, projectName='${projectName}')
dispatch_id = wftk.executeWorkflowFromDB('${workflowName}', scheduler='local')
`;
};
