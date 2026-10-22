import { ApiCall } from './ApiCall';
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

// Reply from polling a run. status is one of running / done / error / not_found.
// error is filled in once the run fails. chunks holds the per-task output segments,
// growing live while the run runs and complete once it is done.
export type PollWorkflowResult = {
  status: string,
  error: string,
  chunks?: WorkflowChunk[] | null,
};

// Every call to /start_workflow and /workflow_status/*. Each method answers
// data, or null when the call failed.
export class WorkflowApi extends ApiCall {
  // Starts a saved workflow run. Returns a token to poll with, or a busy status
  // when the server is already running a workflow.
  async start(projectName: string, doc: ProjectDocument): Promise<StartWorkflowResult | null> {
    return await this.postJson('/start_workflow', { projectName, doc }, { label: 'run workflow', silent: true });
  }

  // Polls a run's status by token. chunks holds the per-task output, live and final.
  async poll(token: string): Promise<PollWorkflowResult | null> {
    return await this.getJson(`/workflow_status/${token}`, {}, { label: 'poll workflow', silent: true });
  }
}
