import { BASEURL } from '../shared/baseurl';

// Reply from starting a run: a token to poll with, or status "busy" when a run
// is already in progress on the server.
export type StartWorkflowResult = {
  token?: string,
  status?: string,
};

// One output segment: the task name (or "__preamble__" / "__between__") and the
// console output captured while that segment was the running one.
export type WorkflowChunk = {
  name: string,
  text: string,
};

// Reply from polling a run. status is one of running / done / error / not_found.
// output and error are filled in once the run is done / failed. chunks holds the
// per-task segments, present only once the run is done (in-process runs).
export type PollWorkflowResult = {
  status: string,
  output: string,
  error: string,
  chunks?: WorkflowChunk[] | null,
};

// Starts a saved workflow run on the server. Returns a token to poll with, or a
// busy status when the server is already running a workflow.
export const startWorkflow = async ({
  projectName,
  workflowName,
}: {
  projectName: string,
  workflowName: string,
}): Promise<StartWorkflowResult> => {
  const response = await fetch(`${BASEURL}/start_workflow`, {
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

// Polls a run's status by token. Once done, output holds the full console output.
export const pollWorkflow = async (token: string): Promise<PollWorkflowResult> => {
  const response = await fetch(`${BASEURL}/workflow_status/${token}`);
  const text = await response.text();
  if (!response.ok) {
    const problem = JSON.parse(text);
    throw new Error(problem.error ?? text);
  }
  return JSON.parse(text);
};
