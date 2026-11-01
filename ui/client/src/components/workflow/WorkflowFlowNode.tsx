import { Box, Stack, Theme, useTheme } from '@mui/material';
import { cyan } from '@mui/material/colors';
import { NodeProps, NodeResizer, useUpdateNodeInternals } from '@xyflow/react';
import { ReactNode, useEffect } from 'react';
import { WorkflowNode } from '../../shared/types';
import { normalizeRequires } from '../../shared/workflow';
import { NodeRunStatus } from './nodeRunStatus';
import { NodeCatalogEntry, nodeTypeIssue, paramsFieldDef } from './nodeCatalog';
import { nodeInputHandleId, nodeOutputHandleId } from './workflowDataflow';
import { WorkflowNodeEditor } from './WorkflowNodeEditor';
import { WorkflowNodeSummary } from './WorkflowNodeSummary';
import { WorkflowNodeDragGrip } from './WorkflowNodeDragGrip';
import { RequiresSide, WorkflowRequiresHandle } from './WorkflowRequiresHandle';

export interface WorkflowFlowNodeData {
  name: string;
  node: WorkflowNode;
  catalog: NodeCatalogEntry[];
  onRename: (newName: string) => void;
  onChange: (node: WorkflowNode) => void;
  onFieldContextMenu: (param: string, x: number, y: number, caret?: number) => void;
  onFieldInlineEdit: (param: string, value: string, caret: number | null, el: HTMLInputElement) => void;
  // How the node did in the last run of this workflow.
  runStatus?: NodeRunStatus;
  // True while this node shows its editor rather than its summary: the pointer
  // is on it, or it is pinned open.
  expanded?: boolean;
  // The icons for the editor's top-right corner, drawn in the order given (the
  // pin, the delete X, ...). The node places them and knows nothing else.
  actionButtons?: ReactNode;
  // True when the open editor should take the node's own space, so the canvas
  // lays the other nodes out around it. Otherwise it is laid over them.
  takeSpace?: boolean;
  // True while the pointer is on this node, which gets a cyan silhouette.
  hovered?: boolean;
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
//
// Away from the pointer it is just a summary card, so a workflow of many nodes
// stays readable. Hovering lays the editor over that card. The card keeps its
// own small size while the editor is open, so the nodes below it never get
// shoved aside by a passing pointer.
export const WorkflowFlowNode = ({ data, selected }: NodeProps) => {
  const { name, node, catalog, onRename, onChange, onFieldContextMenu, onFieldInlineEdit } = data as WorkflowFlowNodeData;
  const runStatus = (data as WorkflowFlowNodeData).runStatus ?? NodeRunStatus.Pending;
  const expanded = (data as WorkflowFlowNodeData).expanded ?? false;
  const actionButtons = (data as WorkflowFlowNodeData).actionButtons;
  const takeSpace = (data as WorkflowFlowNodeData).takeSpace ?? false;
  const hovered = (data as WorkflowFlowNodeData).hovered ?? false;
  const theme = useTheme();

  const params = node.Execution?.input_parameters ?? {};

  const typeIssue = nodeTypeIssue(node, catalog);
  const paramsDef = paramsFieldDef(node, catalog);

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

  // The editor's handles (one per parameter row) mount without changing the
  // node's own size, so ReactFlow has to be told to measure them again.
  const updateNodeInternals = useUpdateNodeInternals();
  useEffect(() => {
    updateNodeInternals(name);
  }, [expanded, name]);

  // A thin cyan silhouette marks the node the pointer is on. It goes on
  // whichever box is the visible frame: the floating editor when there is one,
  // the node itself otherwise.
  const floating = expanded && !takeSpace;
  let silhouetteSx = {};
  if (hovered) {
    silhouetteSx = { outline: `1px solid ${cyan[400]}`, outlineOffset: '2px' };
  }
  let nodeSilhouetteSx = silhouetteSx;
  let editorSilhouetteSx = {};
  if (floating) {
    nodeSilhouetteSx = {};
    editorSilhouetteSx = silhouetteSx;
  }

  // Where the open editor sits: in the node's own box when it takes space, laid
  // over the summary card (and the canvas) when it is only hovered.
  let editorSx = {};
  if (!takeSpace) {
    editorSx = {
      position: 'absolute',
      top: 0,
      left: 0,
      px: 1,
      py: 0.5,
      borderRadius: 1,
      bgcolor: 'background.paper',
      boxShadow: 6,
      zIndex: 10,
    };
  }

  return (
    <Box
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
        ...nodeSilhouetteSx,
        ...runningSx,
      }}
    >
      {/* Drag handles to resize the node; shown while it's selected. Size is
          view-only (React Flow's store), not saved with the workflow. */}
      <NodeResizer isVisible={selected} minWidth={260} minHeight={80} />
      <WorkflowRequiresHandle id={nodeInputHandleId(name)} side={RequiresSide.DependsOn} />
      {/* Outside the editor box on purpose: a grip inside it would be nodrag too. */}
      <WorkflowNodeDragGrip />
      {!takeSpace && (
        <WorkflowNodeSummary
          name={name}
          type={node.type}
          paramNames={Object.keys(params)}
          paramsDef={paramsDef}
          requires={normalizeRequires(node.requires)}
        />
      )}
      {/* The editor. Laid over the summary card while only hovered, so a passing
          pointer leaves the node's measured size, and the nodes below it, alone.
          A node that takes space holds the editor in its own box instead, and
          the canvas lays the others out around it. */}
      {expanded && (
        <Box
          className="nodrag"
          sx={{
            // As wide as the editor's own content needs, as the node itself was
            // before the summary card; never narrower than the card it covers.
            width: 'max-content',
            minWidth: '100%',
            maxWidth: 560,
            ...editorSx,
            ...editorSilhouetteSx,
          }}
        >
        {/* The icon slot: whatever the canvas passed, in that order. */}
        <Stack direction="row" spacing={0.5} sx={{ position: 'absolute', top: -12, right: -12 }}>
          {actionButtons}
        </Stack>
        <WorkflowNodeEditor
          name={name}
          node={node}
          catalog={catalog}
          onRename={onRename}
          onChange={onChange}
          onFieldContextMenu={onFieldContextMenu}
          onFieldInlineEdit={onFieldInlineEdit}
        />
        </Box>
      )}
      <WorkflowRequiresHandle id={nodeOutputHandleId(name)} side={RequiresSide.RequiredBy} />
    </Box>
  );
};
