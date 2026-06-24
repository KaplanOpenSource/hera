import { Box, InputBase, TextField } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { Handle, NodeProps, Position } from '@xyflow/react';
import { useState } from 'react';
import { WorkflowNode } from '../../shared/types';
import { DetailsViewItem, keyForDetailsViewItem } from '../details/DetailsViewItem';
import { WorkflowNodeDeleteButton } from './WorkflowNodeDeleteButton';

export interface WorkflowFlowNodeData {
  name: string;
  node: WorkflowNode;
  onRename: (newName: string) => void;
  onChange: (node: WorkflowNode) => void;
  onDelete: () => void;
  [key: string]: unknown;
}

// Custom ReactFlow node: edits the node name, type, and input parameters in
// place. Delete on hover.
export const WorkflowFlowNode = ({ data, selected }: NodeProps) => {
  const { name, node, onRename, onChange, onDelete } = data as WorkflowFlowNodeData;
  const [draft, setDraft] = useState(name);
  const [hover, setHover] = useState(false);
  const params = node.Execution?.input_parameters ?? {};

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
      onMouseEnter={() => setHover(true)}
      onMouseLeave={() => setHover(false)}
      sx={{
        position: 'relative',
        px: 1,
        py: 0.5,
        minWidth: 260,
        borderRadius: 1,
        bgcolor: 'background.paper',
        border: '1px solid',
        borderColor: selected ? 'primary.main' : 'divider',
      }}
    >
      <Handle type="target" position={Position.Left} />
      {hover && <WorkflowNodeDeleteButton onDelete={onDelete} />}
      <InputBase
        className="nodrag"
        value={draft}
        onChange={(e) => setDraft(e.target.value)}
        onBlur={commit}
        onKeyDown={(e) => { if (e.key === 'Enter') { (e.target as HTMLInputElement).blur(); } }}
        inputProps={{ style: { padding: 0 }, 'aria-label': 'node name' }}
        sx={{ fontSize: 13, fontWeight: 600 }}
      />
      <Box className="nodrag" sx={{ mt: 1 }}>
        <TextField
          label="type"
          size="small"
          fullWidth
          value={node.type ?? ''}
          onChange={(e) => onChange({ ...node, type: e.target.value })}
        />
        <SimpleTreeView
          defaultExpandedItems={[keyForDetailsViewItem('input_parameters')]}
          sx={{ mt: 1 }}
        >
          <DetailsViewItem
            itemKey="input_parameters"
            itemValue={params}
            parentKey={undefined}
            setItemValue={(newVal) => onChange({
              ...node,
              Execution: { ...node.Execution, input_parameters: newVal },
            })}
          />
        </SimpleTreeView>
      </Box>
      <Handle type="source" position={Position.Right} />
    </Box>
  );
};
