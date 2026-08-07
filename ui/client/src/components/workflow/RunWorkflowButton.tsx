import { PlayArrow } from '@mui/icons-material';
import { Button, Dialog, DialogActions, DialogContent, DialogTitle, SxProps, Theme } from '@mui/material';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { runWorkflow } from '../../io/runWorkflow';
import { dismiss, pushError, pushInfo, pushRunning } from '../../io/snackbar';
import { WorkflowLogView } from './log/WorkflowLogView';

// Runs a saved workflow via the server's /run_workflow endpoint. The run is
// synchronous — the button stays busy until the workflow finishes. The captured
// console output is shown in a dialog; errors are surfaced via a snackbar.
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
  const [output, setOutput] = useState<string | null>(null);

  const doRun = async () => {
    const key = pushRunning('run workflow');
    try {
      const { output } = await runWorkflow({ projectName, workflowName });
      setOutput(output ?? '');
      pushInfo(`Workflow "${workflowName}" finished`);
    } catch (e: any) {
      pushError(`run workflow: ${e?.message ?? e}`);
    } finally {
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
      <Dialog open={!!output} onClose={() => setOutput(null)} maxWidth="md" fullWidth>
        <DialogTitle>Workflow "{workflowName}" output</DialogTitle>
        <DialogContent dividers>
          <WorkflowLogView output={output ?? ''} />
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setOutput(null)}>Close</Button>
        </DialogActions>
      </Dialog>
    </>
  );
};
