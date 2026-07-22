import { Chip, Stack, useTheme } from '@mui/material';
import { Handle, Position } from '@xyflow/react';
import { outputHandleId } from './workflowDataflow';

// A single node output: a chip with the output name and a source Handle on the
// node's right edge — the anchor a dataflow line leaves from (id = the output
// name). The row is a fixed height so it lines up with an input parameter row
// (whose tree rows don't margin-collapse the way a plain stack would).
export const WorkflowNodeOutputChip = ({
  nodeName,
  name,
}: {
  nodeName: string,
  name: string,
}) => {
  const theme = useTheme();
  return (
    <Stack direction="row" spacing={0.75} sx={{ alignItems: 'center', justifyContent: 'flex-end', height: '69px' }}>
      <Chip label={name} size="small" variant="outlined" />
      <Handle
        type="source"
        id={outputHandleId(nodeName, name)}
        position={Position.Right}
        style={{ position: 'relative', top: 'auto', right: 'auto', transform: 'none', width: 8, height: 8, background: theme.palette.primary.main, border: 'none' }}
      />
    </Stack>
  );
};
