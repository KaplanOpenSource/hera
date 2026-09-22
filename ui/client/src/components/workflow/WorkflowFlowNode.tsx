import { Autocomplete, Box, InputBase, Stack, TextField, Theme, Typography, useTheme } from '@mui/material';
import { Handle, NodeProps, NodeResizer, Position } from '@xyflow/react';
import { useState } from 'react';
import { WorkflowNode } from '../../shared/types';
import { keyForDetailsViewItem } from '../details/DetailsViewItem';
import { NodeRunStatus } from './nodeRunStatus';
import { NodeCatalogEntry, nodeOutputNames, nodeTypeGroup, nodeTypeIssue, paramsFieldDef } from './nodeCatalog';
import { paramsOnTypeChange } from './nodeTypeParams';
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
  // How the node did in the last run of this workflow.
  runStatus?: NodeRunStatus;
  [key: string]: unknown;
}

// Outline color per run state. Pending has none, so the node keeps its normal border.
const runStatusColor = (status: NodeRunStatus, theme: Theme): string => {
  if (status === NodeRunStatus.Success) {
    return theme.palette.success.main;
  }
  if (status === NodeRunStatus.Failure) {
    return theme.palette.error.main;
  }
  if (status === NodeRunStatus.Running) {
    return theme.palette.info.main;
  }
  return '';
};

// A dashed outline that keeps moving, drawn just outside the node's border.
// Four gradient strips (top, bottom, left, right) slid by one dash each cycle.
const marchingAntsSx = (color: string) => {
  return {
    '@keyframes workflowNodeAnts': {
      to: { backgroundPosition: '26px 0, -26px 100%, 0 -26px, 100% 26px' },
    },
    '&::before': {
      content: '""',
      position: 'absolute',
      inset: '-14px',
      borderRadius: 'inherit',
      pointerEvents: 'none',
      backgroundImage: `linear-gradient(90deg, ${color} 50%, transparent 0), linear-gradient(90deg, ${color} 50%, transparent 0), linear-gradient(0deg, ${color} 50%, transparent 0), linear-gradient(0deg, ${color} 50%, transparent 0)`,
      backgroundSize: '26px 10px, 26px 10px, 10px 26px, 10px 26px',
      backgroundRepeat: 'repeat-x, repeat-x, repeat-y, repeat-y',
      backgroundPosition: '0 0, 0 100%, 0 0, 100% 0',
      animation: 'workflowNodeAnts 0.6s linear infinite',
    },
  };
};

// Custom ReactFlow node: edits the node name, type, and input parameters in
// place. Delete on hover.
export const WorkflowFlowNode = ({ data, selected }: NodeProps) => {
  const { name, node, catalog, onRename, onChange, onDelete, onFieldContextMenu, onFieldInlineEdit } = data as WorkflowFlowNodeData;
  const runStatus = (data as WorkflowFlowNodeData).runStatus ?? NodeRunStatus.Pending;
  const theme = useTheme();
  const [draft, setDraft] = useState(name);
  const [hover, setHover] = useState(false);

  const params = node.Execution?.input_parameters ?? {};

  // Expanded two levels by default: input_parameters and each parameter under it,
  // so nested parameter values are visible without a click.
  const [expandedItems, setExpandedItems] = useState<string[]>(() => {
    const inputsKey = keyForDetailsViewItem(INPUT_PARAMETERS_KEY);
    return [inputsKey, ...Object.keys(params).map(param => keyForDetailsViewItem(param, inputsKey))];
  });

  const typeOptions = catalog.map(entry => entry.type);
  const typeIssue = nodeTypeIssue(node, catalog);
  const paramsDef = paramsFieldDef(node, catalog);
  const outputs = nodeOutputNames(node, catalog);
  const inputsExpanded = expandedItems.includes(keyForDetailsViewItem(INPUT_PARAMETERS_KEY));

  // The run state wins over the type warning and the selection, so a failed or
  // finished node is visible at a glance.
  const statusColor = runStatusColor(runStatus, theme);
  let borderColor = 'divider';
  if (selected) {
    borderColor = 'primary.main';
  }
  if (typeIssue) {
    borderColor = 'warning.main';
  }
  if (statusColor) {
    borderColor = statusColor;
  }
  // Much thicker than the 1px selection border, so the two never look alike.
  let borderWidth = 1;
  if (statusColor) {
    borderWidth = 10;
  }
  let runningSx = {};
  if (runStatus === NodeRunStatus.Running) {
    runningSx = marchingAntsSx(statusColor);
  }

  // Free-form typing keeps the type as-is (custom types stay allowed); picking a
  // known type also seeds its parameters from the catalog.
  const setType = (type: string) => onChange({ ...node, type });
  const pickType = (type: string) => {
    const entry = catalog.find(e => e.type === type);
    onChange({
      ...node,
      type,
      Execution: { ...node.Execution, input_parameters: paramsOnTypeChange(params, entry) },
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
      data-run-status={runStatus}
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
        borderWidth,
        borderColor,
        ...runningSx,
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
