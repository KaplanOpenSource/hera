import { Box, Button, CircularProgress, Dialog, DialogActions, DialogContent, DialogTitle, Typography } from '@mui/material';
import { chunksToText, WorkflowChunk } from '../../../io/runWorkflow';
import { useLogFilterStore } from '../../../stores/useLogFilterStore';
import { chunkedMetrics, flatMetrics, LogMetrics } from './logMetrics';
import { LogToolbar } from './LogToolbar';
import { WorkflowChunkedLog } from './WorkflowChunkedLog';
import { WorkflowLogView } from './WorkflowLogView';

// Shows a workflow run's output, always from the per-task chunks. While the run is
// in progress it shows the flat log (the chunks joined) as it grows with a small
// "running" hint; on finish it shows the grouped per-task view (or an error). The
// log-level filter + copy-all toolbar sits in the bottom action bar, outside the
// scrolling log, so it stays put while the log scrolls.
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
  const visible = useLogFilterStore((state) => { return state.visible; });
  const toggle = useLogFilterStore((state) => { return state.toggle; });

  const chunkList = chunks ?? [];
  const showChunked = !running && chunkList.length > 0;
  // The log area shows while running and after a successful run; on failure it is
  // replaced by the error message, so the toolbar has nothing to act on.
  const showLog = running || !error;

  let metrics: LogMetrics;
  if (showChunked) {
    metrics = chunkedMetrics(chunkList);
  } else {
    metrics = flatMetrics(chunksToText(chunkList));
  }

  return (
    <Dialog open={open} onClose={onClose} maxWidth="md" fullWidth>
      <DialogTitle>Workflow "{workflowName}" output</DialogTitle>
      <DialogContent dividers>
        {!running && error && (
          <Typography color="error" sx={{ whiteSpace: 'pre-wrap', wordBreak: 'break-word' }}>
            {error}
          </Typography>
        )}
        {showLog && (
          <>
            {running && (
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, pb: 1 }}>
                <CircularProgress size={16} />
                <Typography variant="body2" color="text.secondary">Running…</Typography>
              </Box>
            )}
            {showChunked
              ? <WorkflowChunkedLog chunks={chunkList} />
              : <WorkflowLogView output={metrics.fullText} />}
          </>
        )}
      </DialogContent>
      <DialogActions>
        {showLog && (
          <Box sx={{ flexGrow: 1, display: 'flex' }}>
            <LogToolbar counts={metrics.counts} visible={visible} onToggle={toggle} fullText={metrics.fullText} />
          </Box>
        )}
        <Button onClick={onClose}>Close</Button>
      </DialogActions>
    </Dialog>
  );
};
