import { create } from 'zustand';
import { WorkflowChunk } from '../io/runWorkflow';

export enum WorkflowRunStatus {
  Running = 'running',
  Done = 'done',
  Error = 'error',
  Stopped = 'stopped',
}

export type WorkflowRun = {
  token: string,
  status: WorkflowRunStatus,
  error: string,
  // Per-task output segments; the single source of truth, growing live and final.
  chunks: WorkflowChunk[],
};

// A finished run read back from the run history; id names it on the server.
export type OldWorkflowRun = WorkflowRun & {
  id: string,
};

type WorkflowRunStore = {
  // Keyed by workflow name so every button for a workflow shares one run.
  runs: { [workflowName: string]: WorkflowRun },
  // The old run picked in the output tab, shown instead of the latest run.
  oldRuns: { [workflowName: string]: OldWorkflowRun },
  showOldRun: (workflowName: string, run: OldWorkflowRun) => void,
  // Drops the picked old run, so the latest run shows again.
  showLatestRun: (workflowName: string) => void,
  // Called when a run is started; marks the workflow as running.
  startRun: (workflowName: string, token: string) => void,
  // Called by the poller while the run is still going, to show the growing output.
  // Keeps status = running; only updates the chunks.
  setRunChunks: (workflowName: string, chunks: WorkflowChunk[]) => void,
  // Called by the poller when the run ends (done, error or stopped).
  setRunResult: (
    workflowName: string,
    result: { status: WorkflowRunStatus, error: string, chunks: WorkflowChunk[] },
  ) => void,
};

// The run the output tab and the canvas show: the picked old run, else the latest.
export const shownRun = (state: WorkflowRunStore, workflowName: string): WorkflowRun | undefined => {
  return state.oldRuns[workflowName] ?? state.runs[workflowName];
};

// Run state is deliberately not persisted: a token from a previous session would
// be unknown to a freshly started server anyway.
export const useWorkflowRunStore = create<WorkflowRunStore>((set) => {
  return {
    runs: {},
    oldRuns: {},
    showOldRun: (workflowName, run) => {
      return set((state) => {
        return { oldRuns: { ...state.oldRuns, [workflowName]: run } };
      });
    },
    showLatestRun: (workflowName) => {
      return set((state) => {
        const oldRuns = { ...state.oldRuns };
        delete oldRuns[workflowName];
        return { oldRuns };
      });
    },
    startRun: (workflowName, token) => {
      return set((state) => {
        // A new run is shown live, so a picked old run is dropped.
        const oldRuns = { ...state.oldRuns };
        delete oldRuns[workflowName];
        return {
          runs: {
            ...state.runs,
            [workflowName]: { token, status: WorkflowRunStatus.Running, error: '', chunks: [] },
          },
          oldRuns,
        };
      });
    },
    setRunChunks: (workflowName, chunks) => {
      return set((state) => {
        const existing = state.runs[workflowName];
        if (!existing) {
          console.error(`setRunChunks: no run in progress for workflow "${workflowName}"`);
          return {};
        }
        return {
          runs: { ...state.runs, [workflowName]: { ...existing, chunks } },
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
