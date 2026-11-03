import { ApiCall } from './ApiCall';
import { ProjectDocument } from '../shared/types';
import { WorkflowServerStatus } from './WorkflowServerStatus';

// Reply from starting a run: a token to poll with, or status "busy" when a run
// is already in progress on the server.
export type StartWorkflowResult = {
  token?: string,
  status?: WorkflowServerStatus,
};

// One output segment: the task name (or "__preamble__" / "__between__") and the
// console output captured while that segment was the running one.
export type WorkflowChunk = {
  name: string,
  text: string,
};

// Reply from polling a run. status is one of running / done / error / stopped / not_found.
// error is filled in once the run fails. chunks holds the per-task output segments,
// growing live while the run runs and complete once it is done.
export type PollWorkflowResult = {
  status: WorkflowServerStatus,
  error: string,
  chunks?: WorkflowChunk[] | null,
};

// One saved run of a workflow, as the run list shows it. id names it on the server.
export type SavedRunSummary = {
  id: string,
  token: string,
  startTime: string,
  endTime: string,
  status: WorkflowServerStatus,
};

// One saved run with its output.
export type SavedRun = SavedRunSummary & {
  error: string,
  chunks: WorkflowChunk[],
};

// Every call to /start_workflow, /workflow_status/*, /stop_workflow/*, /workflow_runs and /workflow_run. Each method answers
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

  // Stops a running run by token. The poller then sees it as stopped.
  async stop(token: string): Promise<PollWorkflowResult | null> {
    return await this.postJson(`/stop_workflow/${token}`, {}, { label: 'stop workflow', silent: true });
  }

  // The saved runs of a workflow, newest first, without their output.
  async runs(projectName: string, workflowName: string): Promise<SavedRunSummary[] | null> {
    const result = await this.getJson('/workflow_runs', { projectName, workflowName }, { label: 'list workflow runs', silent: true });
    if (!result) {
      return null;
    }
    return result.runs;
  }

  // One saved run of a workflow, with its output.
  async run(projectName: string, workflowName: string, runId: string): Promise<SavedRun | null> {
    return await this.getJson('/workflow_run', { projectName, workflowName, runId }, { label: 'read workflow run' });
  }
}
