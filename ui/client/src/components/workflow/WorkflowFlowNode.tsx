import { Autocomplete, Box, InputBase, Stack, TextField, Typography } from '@mui/material';
import { Handle, NodeProps, Position } from '@xyflow/react';
import { useState } from 'react';
import { WorkflowNode } from '../../shared/types';
import { keyForDetailsViewItem } from '../details/DetailsViewItem';
import { NodeCatalogEntry, nodeOutputNames, nodeTypeGroup, nodeTypeIssue, paramsFieldDef, prefilledParameters } from './nodeCatalog';
import { WorkflowNodeDeleteButton } from './WorkflowNodeDeleteButton';
import { WorkflowNodeInputs } from './WorkflowNodeInputs';
import { WorkflowNodeOutputs } from './WorkflowNodeOutputs';

export interface WorkflowFlowNodeData {
  name: string;
  node: WorkflowNode;
  catalog: NodeCatalogEntry[];
  onRename: (newName: string) => void;
  onChange: (node: WorkflowNode) => void;
  onDelete: () => void;
  onFieldContextMenu: (param: string, x: number, y: number) => void;
  [key: string]: unknown;
}

// Custom ReactFlow node: edits the node name, type, and input parameters in
// place. Delete on hover.
export const WorkflowFlowNode = ({ data, selected }: NodeProps) => {
  const { name, node, catalog, onRename, onChange, onDelete, onFieldContextMenu } = data as WorkflowFlowNodeData;
  const [draft, setDraft] = useState(name);
  const [hover, setHover] = useState(false);
  // Controlled so the input_parameters chevron also shows/hides the outputs.
  const [expandedItems, setExpandedItems] = useState<string[]>([keyForDetailsViewItem('input_parameters')]);
  const params = node.Execution?.input_parameters ?? {};
  const typeOptions = catalog.map(entry => entry.type);
  const typeIssue = nodeTypeIssue(node, catalog);
  const paramsDef = paramsFieldDef(node, catalog);
  const outputs = nodeOutputNames(node, catalog);
  const inputsExpanded = expandedItems.includes(keyForDetailsViewItem('input_parameters'));

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
        borderColor: typeIssue ? 'warning.main' : selected ? 'primary.main' : 'divider',
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
        <Stack direction="row" spacing={1} sx={{ mt: 1, alignItems: 'flex-start' }}>
          <WorkflowNodeInputs
            params={params}
            paramsDef={paramsDef}
            expandedItems={expandedItems}
            onExpandedItemsChange={setExpandedItems}
            onChangeParams={(newVal) => onChange({
              ...node,
              Execution: { ...node.Execution, input_parameters: newVal },
            })}
            onFieldContextMenu={onFieldContextMenu}
          />
          {outputs.length > 0 && (
            <WorkflowNodeOutputs outputs={outputs} expanded={inputsExpanded} />
          )}
        </Stack>
        {typeIssue && (
          <Typography
            className="nodrag"
            variant="caption"
            color="warning.main"
            sx={{ display: 'block', mt: 0.5, userSelect: 'text' }}
          >
            {typeIssue}
          </Typography>
        )}
      </Box>
      <Handle type="source" position={Position.Right} />
    </Box>
  );
};
