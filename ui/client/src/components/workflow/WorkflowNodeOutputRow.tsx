import { Stack, Typography, useTheme } from '@mui/material';
import { Handle, Position } from '@xyflow/react';
import { outputHandleId } from './workflowDataflow';

// A single node output: its name, and a source Handle pushed out to the node's
// right edge — the anchor a dataflow line leaves from (id = the output name).
// An output has this one vertex only, and every output's vertex shares an x.
export const WorkflowNodeOutputRow = ({
  nodeName,
  name,
}: {
  nodeName: string,
  name: string,
}) => {
  const theme = useTheme();
  return (
    <Stack direction="row" spacing={0.75} sx={{ alignItems: 'center', justifyContent: 'space-between' }}>
      <Typography sx={{ whiteSpace: 'nowrap' }}>{name}</Typography>
      <Handle
        type="source"
        id={outputHandleId(nodeName, name)}
        position={Position.Right}
        style={{ position: 'relative', top: 'auto', right: -14, transform: 'none', width: 8, height: 8, background: theme.palette.primary.main, border: 'none' }}
      />
    </Stack>
  );
};
