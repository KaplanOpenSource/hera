import { Box, Button, CircularProgress, Dialog, DialogActions, DialogContent, DialogTitle, Typography } from '@mui/material';
import { WorkflowLogView } from './WorkflowLogView';

// Shows a workflow run's output. Opens immediately with a spinner while the run is
// in progress, then swaps in the classified log (or an error message).
export const WorkflowOutputDialog = ({
  open,
  running,
  output,
  error,
  workflowName,
  onClose,
}: {
  open: boolean,
  running: boolean,
  output: string | null,
  error: string | null,
  workflowName: string,
  onClose: () => void,
}) => {
  return (
    <Dialog open={open} onClose={onClose} maxWidth="md" fullWidth>
      <DialogTitle>Workflow "{workflowName}" output</DialogTitle>
      <DialogContent dividers>
        {running && (
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, py: 3, justifyContent: 'center' }}>
            <CircularProgress size={24} />
            <Typography>Running…</Typography>
          </Box>
        )}
        {!running && error && (
          <Typography color="error" sx={{ whiteSpace: 'pre-wrap', wordBreak: 'break-word' }}>
            {error}
          </Typography>
        )}
        {!running && !error && <WorkflowLogView output={output ?? ''} />}
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose}>Close</Button>
      </DialogActions>
    </Dialog>
  );
};
