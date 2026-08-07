import { PlayArrow } from '@mui/icons-material';
import { SxProps, Theme } from '@mui/material';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { runWorkflow } from '../../io/runWorkflow';
import { dismiss, pushError, pushInfo, pushRunning } from '../../io/snackbar';
import { WorkflowOutputDialog } from './log/WorkflowOutputDialog';

// Runs a saved workflow via the server's /run_workflow endpoint. The run is
// synchronous — the output dialog opens right away and shows a spinner until the
// captured console output (or an error) arrives.
export const RunWorkflowButton = ({
  projectName,
  workflowName,
  disabled,
  disabledReason,
  sx,
}: {
  projectName: string,
  workflowName: string,
  disabled?: boolean,
  disabledReason?: string,
  sx?: SxProps<Theme>,
}) => {
  const [open, setOpen] = useState(false);
  const [running, setRunning] = useState(false);
  const [output, setOutput] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  const doRun = async () => {
    setOpen(true);
    setRunning(true);
    setOutput(null);
    setError(null);
    const key = pushRunning('run workflow');
    try {
      const { output } = await runWorkflow({ projectName, workflowName });
      setOutput(output ?? '');
      pushInfo(`Workflow "${workflowName}" finished`);
    } catch (e: any) {
      setError(e?.message ?? String(e));
      pushError(`run workflow: ${e?.message ?? e}`);
    } finally {
      setRunning(false);
      dismiss(key);
    }
  };

  const title = disabled && disabledReason ? disabledReason : 'Run workflow';
  return (
    <>
      <ButtonTooltip
        title={title}
        aria-label={title}
        disabled={disabled}
        onClick={doRun}
        sx={sx}
      >
        <PlayArrow />
      </ButtonTooltip>
      <WorkflowOutputDialog
        open={open}
        running={running}
        output={output}
        error={error}
        workflowName={workflowName}
        onClose={() => { return setOpen(false); }}
      />
    </>
  );
};
