import { Box, Button, CircularProgress, Dialog, DialogActions, DialogContent, DialogTitle, Typography } from '@mui/material';
import { WorkflowLogView } from './WorkflowLogView';

// Shows a workflow run's output. While the run is in progress it shows the log as
// it grows with a small "running" hint; on finish it shows the final classified log
// (or an error message).
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
        {!running && error && (
          <Typography color="error" sx={{ whiteSpace: 'pre-wrap', wordBreak: 'break-word' }}>
            {error}
          </Typography>
        )}
        {(running || !error) && (
          <>
            {running && (
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, pb: 1 }}>
                <CircularProgress size={16} />
                <Typography variant="body2" color="text.secondary">Running…</Typography>
              </Box>
            )}
            <WorkflowLogView output={output ?? ''} />
          </>
        )}
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose}>Close</Button>
      </DialogActions>
    </Dialog>
  );
};
