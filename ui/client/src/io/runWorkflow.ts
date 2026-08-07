import { BASEURL } from '../shared/baseurl';

export type RunWorkflowResult = {
  dispatch_id: string,
  output: string,
};

// Builds and runs a saved workflow via the server's /run_workflow endpoint,
// returning the dispatch id and the captured console output.
export const runWorkflow = async ({
  projectName,
  workflowName,
}: {
  projectName: string,
  workflowName: string,
}): Promise<RunWorkflowResult> => {
  const response = await fetch(`${BASEURL}/run_workflow`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ projectName, workflowName }),
  });
  const text = await response.text();
  if (!response.ok) {
    const problem = JSON.parse(text);
    throw new Error(problem.error ?? text);
  }
  return JSON.parse(text);
};
