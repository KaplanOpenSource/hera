import { Autocomplete, Box, InputBase, Stack, TextField, Typography } from '@mui/material';
import { Handle, NodeProps, NodeResizer, Position } from '@xyflow/react';
import { useState } from 'react';
import { WorkflowNode } from '../../shared/types';
import { keyForDetailsViewItem } from '../details/DetailsViewItem';
import { NodeCatalogEntry, nodeOutputNames, nodeTypeGroup, nodeTypeIssue, paramsFieldDef, prefilledParameters } from './nodeCatalog';
import { WorkflowNodeDeleteButton } from './WorkflowNodeDeleteButton';
import { nodeInputHandleId, nodeOutputHandleId } from './workflowDataflow';
import { INPUT_PARAMETERS_KEY, WorkflowNodeInputs } from './WorkflowNodeInputs';
import { WorkflowNodeOutputs } from './WorkflowNodeOutputs';

export interface WorkflowFlowNodeData {
  name: string;
  node: WorkflowNode;
  catalog: NodeCatalogEntry[];
  onRename: (newName: string) => void;
  onChange: (node: WorkflowNode) => void;
  onDelete: () => void;
  onFieldContextMenu: (param: string, x: number, y: number, caret?: number) => void;
  onFieldInlineEdit: (param: string, value: string, caret: number | null, el: HTMLInputElement) => void;
  [key: string]: unknown;
}

// Custom ReactFlow node: edits the node name, type, and input parameters in
// place. Delete on hover.
export const WorkflowFlowNode = ({ data, selected }: NodeProps) => {
  const { name, node, catalog, onRename, onChange, onDelete, onFieldContextMenu, onFieldInlineEdit } = data as WorkflowFlowNodeData;
  const [draft, setDraft] = useState(name);
  const [hover, setHover] = useState(false);
  // Controlled so the input_parameters chevron also shows/hides the outputs.
  const [expandedItems, setExpandedItems] = useState<string[]>([keyForDetailsViewItem(INPUT_PARAMETERS_KEY)]);
  const params = node.Execution?.input_parameters ?? {};
  const typeOptions = catalog.map(entry => entry.type);
  const typeIssue = nodeTypeIssue(node, catalog);
  const paramsDef = paramsFieldDef(node, catalog);
  const outputs = nodeOutputNames(node, catalog);
  const inputsExpanded = expandedItems.includes(keyForDetailsViewItem(INPUT_PARAMETERS_KEY));

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
        // Fill the node wrapper so a drag-resize (which sizes the wrapper) grows
        // this box with it; before any resize the wrapper shrink-wraps content.
        width: '100%',
        height: '100%',
        borderRadius: 1,
        bgcolor: 'background.paper',
        border: '1px solid',
        borderColor: typeIssue ? 'warning.main' : selected ? 'primary.main' : 'divider',
      }}
    >
      {/* Drag handles to resize the node; shown while it's selected. Size is
          view-only (React Flow's store), not saved with the workflow. */}
      <NodeResizer isVisible={selected} minWidth={260} minHeight={80} />
      <Handle type="target" id={nodeInputHandleId(name)} position={Position.Left} />
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
            nodeName={name}
            params={params}
            paramsDef={paramsDef}
            expandedItems={expandedItems}
            onExpandedItemsChange={setExpandedItems}
            onChangeParams={(newVal) => onChange({
              ...node,
              Execution: { ...node.Execution, input_parameters: newVal },
            })}
            onFieldContextMenu={onFieldContextMenu}
            onFieldInlineEdit={onFieldInlineEdit}
          />
          {outputs.length > 0 && (
            <WorkflowNodeOutputs nodeName={name} outputs={outputs} expanded={inputsExpanded} />
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
      <Handle type="source" id={nodeOutputHandleId(name)} position={Position.Right} />
    </Box>
  );
};
