import { useEffect, useRef } from 'react';
import { SnackbarKey } from 'notistack';
import { pollWorkflow } from '../../io/runWorkflow';
import { dismiss, pushError, pushInfo, pushRunning } from '../../io/snackbar';
import { useWorkflowRunStore, WorkflowRunStatus } from '../../stores/useWorkflowRunStore';

// How often to poll a running workflow for its status.
const POLL_MS = 500;

// Mounted once (in App). Watches the workflow-run store for a running run and
// polls it until it finishes, writing the result back to the store so every run
// button for that workflow sees it. The polling useEffect lives here, out of the
// buttons, so the run survives a button unmounting. The server runs one workflow
// at a time, so there is at most one running run to poll.
export const WorkflowRunPoller = () => {
  const runs = useWorkflowRunStore((state) => { return state.runs; });
  const setRunResult = useWorkflowRunStore((state) => { return state.setRunResult; });
  const setRunOutput = useWorkflowRunStore((state) => { return state.setRunOutput; });

  const runningEntry = Object.entries(runs).find(([, run]) => {
    return run.status === WorkflowRunStatus.Running;
  });
  const workflowName = runningEntry ? runningEntry[0] : null;
  const token = runningEntry ? runningEntry[1].token : null;

  const snackbarKeyRef = useRef<SnackbarKey | null>(null);

  useEffect(() => {
    if (!token || !workflowName) {
      return;
    }
    let cancelled = false;
    let timer: ReturnType<typeof setTimeout>;
    snackbarKeyRef.current = pushRunning('run workflow');

    const clearRunningSnackbar = () => {
      if (snackbarKeyRef.current) {
        dismiss(snackbarKeyRef.current);
        snackbarKeyRef.current = null;
      }
    };

    const poll = async () => {
      let result;
      try {
        result = await pollWorkflow(token);
      } catch (e: any) {
        if (cancelled) {
          return;
        }
        setRunResult(workflowName, { status: WorkflowRunStatus.Error, output: '', error: e?.message ?? String(e) });
        pushError(`run workflow: ${e?.message ?? e}`);
        clearRunningSnackbar();
        return;
      }
      if (cancelled) {
        return;
      }
      if (result.status === 'running') {
        // Show the output as it grows, without ending the run.
        setRunOutput(workflowName, result.output ?? '');
        timer = setTimeout(poll, POLL_MS);
      } else if (result.status === 'done') {
        setRunResult(workflowName, { status: WorkflowRunStatus.Done, output: result.output ?? '', error: '' });
        pushInfo(`Workflow "${workflowName}" finished`);
        clearRunningSnackbar();
      } else if (result.status === 'error') {
        const message = result.error || 'Workflow failed';
        setRunResult(workflowName, { status: WorkflowRunStatus.Error, output: '', error: message });
        pushError(`run workflow: ${message}`);
        clearRunningSnackbar();
      } else {
        // not_found: the server restarted or a newer run overwrote the slot.
        const message = 'The run was lost (the server may have restarted).';
        setRunResult(workflowName, { status: WorkflowRunStatus.Error, output: '', error: message });
        pushError(`run workflow: ${message}`);
        clearRunningSnackbar();
      }
    };
    poll();

    return () => {
      cancelled = true;
      clearTimeout(timer);
      clearRunningSnackbar();
    };
  }, [token, workflowName, setRunResult, setRunOutput]);

  return null;
};
