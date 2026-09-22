import { Box, useTheme } from '@mui/material';
import { Background, Controls, Edge, ReactFlow, ReactFlowProvider } from '@xyflow/react';
import '@xyflow/react/dist/style.css';
import { ReactNode, useMemo, useState } from 'react';
import { WorkflowBlock, WorkflowNode } from '../../shared/types';
import { NodeCatalogEntry } from './nodeCatalog';
import { WorkflowCanvasEdits } from './WorkflowCanvasEdits';
import { WorkflowCanvasToolbar } from './WorkflowCanvasToolbar';
import { WorkflowDisplayEdges } from './WorkflowDisplayEdges';
import { WorkflowReferences } from './WorkflowReferences';
import { FIT_MAX_ZOOM } from './WorkflowViewport';
import { WorkflowContextMenu, WorkflowContextMenuKind, WorkflowContextMenuTarget } from './WorkflowContextMenu';
import { WorkflowFlowNode } from './WorkflowFlowNode';
import { WorkflowRequiresEdge } from './WorkflowRequiresEdge';
import { buildWorkflowEdges } from './workflowEdges';
import { buildDataflowEdges } from './workflowDataflow';
import { WorkflowInlineReference } from './WorkflowInlineReference';
import { displayFlowNodes, flowStructureKey } from './workflowFlowNodesUtils';
import { useInlineReference } from './useInlineReference';
import { useWorkflowCanvasNodes } from './useWorkflowCanvasNodes';
import { NodeRunStatusMap } from './nodeRunStatus';

// Defined once (module scope) so ReactFlow doesn't warn about changing types.
const NODE_TYPES = { workflow: WorkflowFlowNode };
// Dataflow edges reuse the same removable-edge component (X button at midpoint).
const EDGE_TYPES = { requires: WorkflowRequiresEdge, dataflow: WorkflowRequiresEdge };

interface WorkflowGraphProps {
  catalog: NodeCatalogEntry[];
  nodeNames: string[];
  nodes: { [name: string]: WorkflowNode };
  selectedNode?: string;
  // How each node did in the last run, keyed by node name. Missing means pending.
  nodeStatuses?: NodeRunStatusMap;
  // A request from the output tab to centre on one node.
  focus?: { nodeName: string, seq: number };
  // The node the pointer is over, or undefined when it left one.
  onHoverNode?: (name: string | undefined) => void;
  // Extra controls rendered in the canvas's top-right corner (e.g. a run button).
  actionButtons?: ReactNode;
  onSelectNode: (name: string | undefined) => void;
  onAddNode: () => void;
  onApplyTemplate: (block: WorkflowBlock) => void;
  onRenameNode: (oldName: string, newName: string) => void;
  onSetNode: (name: string, node: WorkflowNode) => void;
  onAddRequire: (source: string, target: string) => void;
  onRemoveRequire: (source: string, target: string) => void;
  onDeleteNode: (name: string) => void;
}

// Node graph view of a workflow: one node per workflow node, edges from the
// `requires` field. Nodes are draggable, their names editable inline, and the
// on-canvas button adds a node. Clicking a node selects it for editing.
const WorkflowGraph = ({
  catalog,
  nodeNames,
  nodes,
  selectedNode,
  nodeStatuses,
  focus,
  onHoverNode,
  actionButtons,
  onSelectNode,
  onAddNode,
  onApplyTemplate,
  onRenameNode,
  onSetNode,
  onAddRequire,
  onRemoveRequire,
  onDeleteNode,
}: WorkflowGraphProps) => {
  const theme = useTheme();
  const [menu, setMenu] = useState<WorkflowContextMenuTarget | null>(null);
  const [hoveredEdge, setHoveredEdge] = useState<string | null>(null);
  // Dataflow dependencies (an input referencing another node's output) — used
  // both to draw the lines and to order the columns, alongside `requires`.
  const dataflowDeps = buildDataflowEdges(nodeNames, nodes, catalog);

  // Which node outputs each field may point at (the reference menus ask this).
  const references = new WorkflowReferences(nodeNames, nodes, catalog);

  // What a canvas gesture (a dragged line, a menu command) does to the workflow.
  const edits = new WorkflowCanvasEdits(nodeNames, nodes, { onSetNode, onAddRequire, onRemoveRequire });

  // Where the nodes sit on the canvas and where the view looks.
  const { rfNodes, onNodesChange, containerRef } = useWorkflowCanvasNodes({ nodeNames, nodes, dataflowDeps, focus });

  // The {…} reference autocomplete on a node's parameter fields.
  const { inline, closeInline, handleInlineEdit, pickInline } = useInlineReference({ nodes, references, onSetNode });

  const displayNodes = displayFlowNodes({
    rfNodes,
    nodes,
    catalog,
    nodeStatuses,
    selectedNode,
    handlers: {
      onRename: onRenameNode,
      onChange: onSetNode,
      onDelete: onDeleteNode,
      onFieldContextMenu: (name, param, x, y, caret) =>
        setMenu({ kind: WorkflowContextMenuKind.Field, node: name, param, x, y, caret }),
      onFieldInlineEdit: (name, param, value, caret, el) => handleInlineEdit(name, param, value, caret, el),
    },
  });

  const rfEdges = useMemo<Edge[]>(() => buildWorkflowEdges(nodeNames, nodes), [flowStructureKey(nodeNames, nodes, dataflowDeps)]);

  const displayEdges = WorkflowDisplayEdges.hovering(hoveredEdge)
    .withRequires(rfEdges, onRemoveRequire)
    .withDataflow(dataflowDeps, theme.palette.primary.main, id => edits.removeDataflowEdge(id))
    .all();

  // While a field's menu is open, the other nodes that produce outputs — the
  // options for its "Reference another node's output…" autocomplete steps.
  let referenceOptions: { node: string, outputs: string[] }[] = [];
  if (menu?.kind === WorkflowContextMenuKind.Field) {
    referenceOptions = references.optionsFor(menu.node);
  }

  return (
    <Box ref={containerRef} sx={{ flex: 1, minHeight: 200, ml: -2, mr: -2, mb: -2, borderTop: '1px solid', borderColor: 'divider' }}>
      <ReactFlow
        colorMode={theme.palette.mode}
        nodes={displayNodes}
        edges={displayEdges}
        nodeTypes={NODE_TYPES}
        edgeTypes={EDGE_TYPES}
        onNodesChange={onNodesChange}
        onConnect={connection => edits.connect(connection)}
        onEdgesDelete={deleted => edits.removeEdges(deleted)}
        isValidConnection={connection => edits.canConnect(connection)}
        fitView
        fitViewOptions={{ maxZoom: FIT_MAX_ZOOM }}
        onNodeClick={(_e, node) => onSelectNode(node.id)}
        onPaneClick={() => { onSelectNode(undefined); closeInline(); }}
        onEdgeMouseEnter={(_e, edge) => setHoveredEdge(edge.id)}
        onEdgeMouseLeave={() => setHoveredEdge(null)}
        onNodeContextMenu={(event, node) => {
          event.preventDefault();
          setMenu({ kind: WorkflowContextMenuKind.Node, name: node.id, x: event.clientX, y: event.clientY });
        }}
        onNodeMouseEnter={(_event, node) => { return onHoverNode?.(node.id); }}
        onNodeMouseLeave={() => { return onHoverNode?.(undefined); }}
        onEdgeContextMenu={(event, edge) => {
          event.preventDefault();
          setMenu({ kind: WorkflowContextMenuKind.Edge, source: edge.source, target: edge.target, x: event.clientX, y: event.clientY });
        }}
      >
        <WorkflowCanvasToolbar
          actionButtons={actionButtons}
          onAddNode={onAddNode}
          onApplyTemplate={onApplyTemplate}
        />
        <Background />
        <Controls />
      </ReactFlow>
      <WorkflowContextMenu
        menu={menu}
        referenceOptions={referenceOptions}
        onClose={() => setMenu(null)}
        onDeleteNode={onDeleteNode}
        onDeleteField={(nodeName, param) => edits.deleteField(nodeName, param)}
        onRemoveRequire={onRemoveRequire}
        onReferenceOutput={(nodeName, param, source, output, caret) => edits.referenceOutput(nodeName, param, source, output, caret)}
      />
      <WorkflowInlineReference
        anchorEl={inline?.anchorEl ?? null}
        options={inline?.options ?? []}
        onPick={pickInline}
        onClose={closeInline}
      />
    </Box>
  );
};

// ReactFlowProvider supplies the store that hooks like useNodesInitialized read,
// so the inner component must live inside it.
export const WorkflowGraphWrapper = (props: WorkflowGraphProps) => {
  return (
    <ReactFlowProvider>
      <WorkflowGraph {...props} />
    </ReactFlowProvider>
  );
};
