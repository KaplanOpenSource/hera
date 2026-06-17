import { Box, InputBase, Typography } from '@mui/material';
import { Handle, NodeProps, Position } from '@xyflow/react';
import { useState } from 'react';

export interface WorkflowFlowNodeData {
  name: string;
  type?: string;
  onRename: (newName: string) => void;
  [key: string]: unknown;
}

// Custom ReactFlow node: shows the node's name (editable inline) and its type.
export const WorkflowFlowNode = ({ data, selected }: NodeProps) => {
  const { name, type, onRename } = data as WorkflowFlowNodeData;
  const [draft, setDraft] = useState(name);

  const commit = () => {
    const next = draft.trim();
    if (next && next !== name) {
      onRename(next);
    } else {
      setDraft(name);
    }
  };

  return (
    <Box
      sx={{
        px: 1,
        py: 0.5,
        minWidth: 120,
        borderRadius: 1,
        bgcolor: 'background.paper',
        border: '1px solid',
        borderColor: selected ? 'primary.main' : 'divider',
      }}
    >
      <Handle type="target" position={Position.Left} />
      <InputBase
        className="nodrag"
        value={draft}
        onChange={(e) => setDraft(e.target.value)}
        onBlur={commit}
        onKeyDown={(e) => { if (e.key === 'Enter') { (e.target as HTMLInputElement).blur(); } }}
        inputProps={{ style: { padding: 0 }, 'aria-label': 'node name' }}
        sx={{ fontSize: 13, fontWeight: 600 }}
      />
      <Typography sx={{ fontSize: 10, opacity: 0.7 }}>{type || '—'}</Typography>
      <Handle type="source" position={Position.Right} />
    </Box>
  );
};
