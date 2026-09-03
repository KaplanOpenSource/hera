import { create } from 'zustand';
import { WorkflowChunk } from '../io/runWorkflow';

export enum WorkflowRunStatus {
  Running = 'running',
  Done = 'done',
  Error = 'error',
}

export type WorkflowRun = {
  token: string,
  status: WorkflowRunStatus,
  output: string,
  error: string,
  // Per-task output segments; set when the run finishes (null while running).
  chunks?: WorkflowChunk[] | null,
};

type WorkflowRunStore = {
  // Keyed by workflow name so every button for a workflow shares one run.
  runs: { [workflowName: string]: WorkflowRun },
  // Called when a run is started; marks the workflow as running.
  startRun: (workflowName: string, token: string) => void,
  // Called by the poller while the run is still going, to show the growing output.
  // Keeps status = running; only updates the output text.
  setRunOutput: (workflowName: string, output: string) => void,
  // Called by the poller when the run finishes (done or error).
  setRunResult: (
    workflowName: string,
    result: { status: WorkflowRunStatus, output: string, error: string, chunks?: WorkflowChunk[] | null },
  ) => void,
};

// Run state is deliberately not persisted: a token from a previous session would
// be unknown to a freshly started server anyway.
export const useWorkflowRunStore = create<WorkflowRunStore>((set) => {
  return {
    runs: {},
    startRun: (workflowName, token) => {
      return set((state) => {
        return {
          runs: {
            ...state.runs,
            [workflowName]: { token, status: WorkflowRunStatus.Running, output: '', error: '' },
          },
        };
      });
    },
    setRunOutput: (workflowName, output) => {
      return set((state) => {
        const existing = state.runs[workflowName];
        if (!existing) {
          console.error(`setRunOutput: no run in progress for workflow "${workflowName}"`);
          return {};
        }
        return {
          runs: { ...state.runs, [workflowName]: { ...existing, output } },
        };
      });
    },
    setRunResult: (workflowName, result) => {
      return set((state) => {
        const existing = state.runs[workflowName];
        if (!existing) {
          // The poller only finishes a run that startRun added, so this should
          // never happen. Surface it (don't swallow) but leave state unchanged
          // rather than crash the poll flow.
          console.error(`setRunResult: no run in progress for workflow "${workflowName}"`);
          return {};
        }
        return {
          runs: { ...state.runs, [workflowName]: { ...existing, ...result } },
        };
      });
    },
  };
});
