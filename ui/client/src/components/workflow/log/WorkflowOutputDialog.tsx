import { Box, Button, CircularProgress, Dialog, DialogActions, DialogContent, DialogTitle, Typography } from '@mui/material';
import { chunksToText, WorkflowChunk } from '../../../io/runWorkflow';
import { WorkflowChunkedLog } from './WorkflowChunkedLog';
import { WorkflowLogView } from './WorkflowLogView';

// Shows a workflow run's output, always from the per-task chunks. While the run is
// in progress it shows the flat log (the chunks joined) as it grows with a small
// "running" hint; on finish it shows the grouped per-task view (or an error).
export const WorkflowOutputDialog = ({
  open,
  running,
  chunks,
  error,
  workflowName,
  onClose,
}: {
  open: boolean,
  running: boolean,
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
              : <WorkflowLogView output={chunksToText(chunks)} />}
          </>
        )}
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose}>Close</Button>
      </DialogActions>
    </Dialog>
  );
};
