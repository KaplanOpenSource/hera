import { Box, Button, CircularProgress, Dialog, DialogActions, DialogContent, DialogTitle, Typography } from '@mui/material';
import { WorkflowChunk } from '../../../io/runWorkflow';
import { WorkflowChunkedLog } from './WorkflowChunkedLog';
import { WorkflowLogView } from './WorkflowLogView';

// Shows a workflow run's output. While the run is in progress it shows the flat log
// as it grows with a small "running" hint; on finish, if per-task segments are
// available it shows the grouped view, otherwise the final flat classified log
// (or an error message).
export const WorkflowOutputDialog = ({
  open,
  running,
  output,
  chunks,
  error,
  workflowName,
  onClose,
}: {
  open: boolean,
  running: boolean,
  output: string | null,
  chunks?: WorkflowChunk[] | null,
  error: string | null,
  workflowName: string,
  onClose: () => void,
}) => {
  const showChunked = !running && chunks != null && chunks.length > 0;
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
            {showChunked
              ? <WorkflowChunkedLog chunks={chunks} />
              : <WorkflowLogView output={output ?? ''} />}
          </>
        )}
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose}>Close</Button>
      </DialogActions>
    </Dialog>
  );
};
