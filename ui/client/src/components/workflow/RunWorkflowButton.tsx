import { PlayArrow } from '@mui/icons-material';
import { Box, Button, Dialog, DialogActions, DialogContent, DialogTitle, SxProps, Theme } from '@mui/material';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { fetchPython } from '../../io/fetchPython';
import { pushInfo } from '../../io/snackbar';
import { buildRunWorkflowCode } from './buildRunWorkflowCode';

// Runs a saved workflow by building and executing it from the DB. The run is
// synchronous — the button stays busy until the workflow finishes (errors are
// surfaced by fetchPython's snackbar). The captured console output is shown in a dialog.
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

  const runWorkflow = async () => {
    const { data } = await fetchPython({
      results: ['dispatch_id', 'output'],
      label: 'run workflow',
      code: buildRunWorkflowCode({ projectName, workflowName }),
    });
    if (!data) {
      return;
    }
    setOutput(data.output ?? '');
    pushInfo(`Workflow "${workflowName}" finished`);
  };

  const title = disabled && disabledReason ? disabledReason : 'Run workflow';
  return (
    <>
      <ButtonTooltip
        title={title}
        aria-label={title}
        disabled={disabled}
        onClick={runWorkflow}
        sx={sx}
      >
        <PlayArrow />
      </ButtonTooltip>
      <Dialog open={!!output} onClose={() => setOutput(null)} maxWidth="md" fullWidth>
        <DialogTitle>Workflow "{workflowName}" output</DialogTitle>
        <DialogContent dividers>
          <Box
            component="pre"
            sx={{ m: 0, whiteSpace: 'pre-wrap', wordBreak: 'break-word', fontFamily: 'monospace', fontSize: 12 }}
          >
            {output}
          </Box>
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setOutput(null)}>Close</Button>
        </DialogActions>
      </Dialog>
    </>
  );
};
