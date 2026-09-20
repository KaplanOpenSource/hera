import { BASEURL } from '../shared/baseurl';
import { ProjectDocument } from '../shared/types';

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

// Joins per-task chunks into one flat text, in run order (for a plain-text view).
export const chunksToText = (
  chunks?: WorkflowChunk[] | null,
): string => {
  return (chunks ?? []).map((chunk) => { return chunk.text; }).join('');
};

// Reply from polling a run. status is one of running / done / error / not_found.
// error is filled in once the run fails. chunks holds the per-task output segments,
// growing live while the run runs and complete once it is done.
export type PollWorkflowResult = {
  status: string,
  error: string,
  chunks?: WorkflowChunk[] | null,
};

// Starts a saved workflow run on the server. Returns a token to poll with, or a
// busy status when the server is already running a workflow.
export const startWorkflow = async ({
  projectName,
  doc,
}: {
  projectName: string,
  // The whole workflow document. Sent so the server builds from it, no DB lookup.
  doc: ProjectDocument,
}): Promise<StartWorkflowResult> => {
  const response = await fetch(`${BASEURL}/start_workflow`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ projectName, doc }),
  });
  const text = await response.text();
  if (!response.ok) {
    const problem = JSON.parse(text);
    throw new Error(problem.error ?? text);
  }
  return JSON.parse(text);
};

// Polls a run's status by token. chunks holds the per-task output, live and final.
export const pollWorkflow = async (token: string): Promise<PollWorkflowResult> => {
  const response = await fetch(`${BASEURL}/workflow_status/${token}`);
  const text = await response.text();
  if (!response.ok) {
    const problem = JSON.parse(text);
    throw new Error(problem.error ?? text);
  }
  return JSON.parse(text);
};
