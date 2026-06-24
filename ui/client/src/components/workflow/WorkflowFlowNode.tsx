import { Autocomplete, Box, InputBase, TextField, Typography } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { Handle, NodeProps, Position } from '@xyflow/react';
import { useState } from 'react';
import { WorkflowNode } from '../../shared/types';
import { DetailsViewItem, keyForDetailsViewItem } from '../details/DetailsViewItem';
import { NodeCatalogEntry, nodeTypeGroup, prefilledParameters, validateNode } from './nodeCatalog';
import { WorkflowNodeDeleteButton } from './WorkflowNodeDeleteButton';

export interface WorkflowFlowNodeData {
  name: string;
  node: WorkflowNode;
  catalog: NodeCatalogEntry[];
  onRename: (newName: string) => void;
  onChange: (node: WorkflowNode) => void;
  onDelete: () => void;
  [key: string]: unknown;
}

// Custom ReactFlow node: edits the node name, type, and input parameters in
// place. Delete on hover.
export const WorkflowFlowNode = ({ data, selected }: NodeProps) => {
  const { name, node, catalog, onRename, onChange, onDelete } = data as WorkflowFlowNodeData;
  const [draft, setDraft] = useState(name);
  const [hover, setHover] = useState(false);
  const params = node.Execution?.input_parameters ?? {};
  const typeOptions = catalog.map(entry => entry.type);
  const problems = validateNode(node, catalog);

  // Free-form typing keeps the type as-is (custom types stay allowed); picking a
  // known type also seeds its parameters from the catalog.
  const setType = (type: string) => onChange({ ...node, type });
  const pickType = (type: string) => {
    const entry = catalog.find(e => e.type === type);
    onChange({
      ...node,
      type,
      Execution: { ...node.Execution, input_parameters: prefilledParameters(entry, params) },
    });
  };

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
        borderColor: problems.length > 0 ? 'warning.main' : selected ? 'primary.main' : 'divider',
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
        <Autocomplete
          className="nodrag"
          freeSolo
          size="small"
          options={typeOptions}
          groupBy={(option) => nodeTypeGroup(option)}
          inputValue={node.type ?? ''}
          onInputChange={(_e, value, reason) => {
            if (reason === 'input') {
              setType(value);
            } else if (reason === 'clear') {
              setType('');
            }
          }}
          onChange={(_e, value) => pickType(typeof value === 'string' ? value : value ?? '')}
          renderInput={(inputParams) => <TextField {...inputParams} label="type" fullWidth />}
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
        {problems.length > 0 && (
          <Typography
            className="nodrag"
            variant="caption"
            color="warning.main"
            sx={{ display: 'block', mt: 0.5, whiteSpace: 'pre-wrap', userSelect: 'text' }}
          >
            {problems.join('\n')}
          </Typography>
        )}
      </Box>
      <Handle type="source" position={Position.Right} />
    </Box>
  );
};
