import { Box } from '@mui/material';
import { WorkflowChunk } from '../../../io/runWorkflow';
import { WorkflowChunkLog } from './WorkflowChunkLog';

// Renders a run's output as a stack of per-task cards (one WorkflowChunkLog per
// chunk). The shared filter + copy-all toolbar lives in the dialog's action bar.
export const WorkflowChunkedLog = ({
  chunks,
}: {
  chunks: WorkflowChunk[],
}) => {
  return (
    <Box sx={{ fontFamily: 'monospace', fontSize: 12 }}>
      {chunks.map((chunk, chunkIndex) => {
        return <WorkflowChunkLog key={chunkIndex} chunk={chunk} />;
      })}
    </Box>
  );
};
