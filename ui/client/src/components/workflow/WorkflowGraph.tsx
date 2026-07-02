import { Add } from '@mui/icons-material';
import { Box, IconButton, Tooltip } from '@mui/material';
import { Background, Connection, Controls, Edge, MarkerType, Node, Panel, ReactFlow, ReactFlowProvider, useNodesInitialized, useNodesState, useReactFlow } from '@xyflow/react';
import '@xyflow/react/dist/style.css';
import { useEffect, useMemo, useRef, useState } from 'react';
import { WorkflowNode } from '../../shared/types';
import { NodeCatalogEntry } from './nodeCatalog';
import { WorkflowContextMenu, WorkflowContextMenuKind, WorkflowContextMenuTarget } from './WorkflowContextMenu';
import { WorkflowFlowNode } from './WorkflowFlowNode';
import { WorkflowRequiresEdge } from './WorkflowRequiresEdge';
import { buildWorkflowEdges, isValidConnection as isValidConnectionPure } from './workflowEdges';
import { buildDataflowEdges, clearInputReference, parseDataflowConnection, parseDataflowEdgeId, setInputReference } from './workflowDataflow';
import { WorkflowLayout } from './WorkflowLayout';
import { computeLayers } from './workflowGeometry';

// Defined once (module scope) so ReactFlow doesn't warn about changing types.
const NODE_TYPES = { workflow: WorkflowFlowNode };
// Dataflow edges reuse the same removable-edge component (X button at midpoint).
const EDGE_TYPES = { requires: WorkflowRequiresEdge, dataflow: WorkflowRequiresEdge };

interface WorkflowGraphProps {
  catalog: NodeCatalogEntry[];
  nodeNames: string[];
  nodes: { [name: string]: WorkflowNode };
  selectedNode?: string;
  onSelectNode: (name: string | undefined) => void;
  onAddNode: () => void;
  onRenameNode: (oldName: string, newName: string) => void;
  onSetNode: (name: string, node: WorkflowNode) => void;
  onAddRequire: (source: string, target: string) => void;
  onRemoveRequire: (source: string, target: string) => void;
  onDeleteNode: (name: string) => void;
}

// Node graph view of a workflow: one node per workflow node, edges from the
// `requires` field. Nodes are draggable, their names editable inline, and the
// on-canvas button adds a node. Clicking a node selects it for editing.
const WorkflowGraphInner = ({
  catalog,
  nodeNames,
  nodes,
  selectedNode,
  onSelectNode,
  onAddNode,
  onRenameNode,
  onSetNode,
  onAddRequire,
  onRemoveRequire,
  onDeleteNode,
}: WorkflowGraphProps) => {
  const [menu, setMenu] = useState<WorkflowContextMenuTarget | null>(null);
  const [hoveredEdge, setHoveredEdge] = useState<string | null>(null);
  const { fitView, getViewport, setViewport, getNode, setCenter } = useReactFlow();
  const nodesInitialized = useNodesInitialized();
  const containerRef = useRef<HTMLDivElement>(null);
  const prevHeightRef = useRef<number | null>(null);
  // What to do once nodes are measured after a structure change: 'all' fits the
  // whole graph (initial load / bulk), a node name pans to focus that newly
  // added node while keeping the current zoom. prevNames detects the change.
  const pendingRef = useRef<'all' | string | null>(null);
  const prevNamesRef = useRef<string[]>([]);
  // useNodesState owns only position and identity; the structure effect rebuilds
  // it (preserving dragged positions) when nodes are added/removed/reordered.
  const [rfNodes, setRfNodes, onNodesChange] = useNodesState<Node>([]);

  // Dataflow dependencies (an input referencing another node's output) — used
  // both to draw the lines and to order the columns, alongside `requires`.
  const dataflowDeps = buildDataflowEdges(nodeNames, nodes, catalog);

  // A signature of the workflow structure (names, types, requires, and dataflow
  // links) so the graph only rebuilds when the structure changes — not on drag.
  const structureKey = JSON.stringify([
    nodeNames.map(name => [name, nodes[name]?.type, nodes[name]?.requires]),
    dataflowDeps.map(edge => [edge.source, edge.target]),
  ]);

  useEffect(() => {
    const layout = WorkflowLayout.stacked(nodeNames, nodes, dataflowDeps).positions();
    setRfNodes(prev => {
      const prevPos = new Map(prev.map(n => [n.id, n.position]));
      return nodeNames.map(name => ({
        id: name,
        type: 'workflow',
        position: prevPos.get(name) ?? layout[name],
        data: {},
      }));
    });
  }, [structureKey]);

  // Re-stack whenever nodes are added or removed (and on mount) — but not on
  // every drag or param edit. A single added node is focused (pan to it, keep
  // zoom) so it isn't lost off-screen when fitting all would zoom out too far;
  // initial load / bulk changes fit the whole graph; removes fit what remains.
  // Re-stack when the column assignment changes: nodes added/removed, or a
  // dependency (requires or a dataflow reference) moves a node to another column.
  const layerKey = JSON.stringify(computeLayers(nodeNames, nodes, dataflowDeps));
  useEffect(() => {
    const layout = WorkflowLayout.stacked(nodeNames, nodes, dataflowDeps).positions();
    setRfNodes(prev => prev.map(node => ({ ...node, position: layout[node.id] ?? node.position })));
    const isInitial = prevNamesRef.current.length === 0;
    const added = nodeNames.filter(n => !prevNamesRef.current.includes(n));
    const removed = prevNamesRef.current.filter(n => !nodeNames.includes(n));
    prevNamesRef.current = nodeNames;
    if (isInitial) {
      pendingRef.current = 'all';
    } else if (added.length > 0 && removed.length === 0) {
      pendingRef.current = added.length === 1 ? added[0] : 'all';
    } else {
      requestAnimationFrame(() => fitView({ duration: 300 }));
    }
  }, [layerKey]);

  // Once nodes are measured after a structure change, fit the whole graph or pan
  // to focus the newly added node (keeping the current zoom).
  useEffect(() => {
    if (!nodesInitialized || pendingRef.current === null) {
      return;
    }
    const pending = pendingRef.current;
    pendingRef.current = null;
    if (pending === 'all') {
      fitView({ duration: 300 });
      return;
    }
    const node = getNode(pending);
    if (node) {
      const x = node.position.x + (node.measured?.width ?? 0) / 2;
      const y = node.position.y + (node.measured?.height ?? 0) / 2;
      setCenter(x, y, { zoom: getViewport().zoom, duration: 300 });
    }
  }, [nodesInitialized]);

  // Once nodes are measured, push down only the ones that overlap within their
  // column (using real measured heights) — so growing a node, e.g. by picking a
  // type with more parameters, shoves the nodes below it instead of overlapping
  // them, while leaving every non-colliding position (including drags) untouched.
  const measuredKey = JSON.stringify(rfNodes.map(node => [node.id, Math.round(node.measured?.height ?? 0)]));
  useEffect(() => {
    const fixed = WorkflowLayout.fromFlowNodes(rfNodes, nodeNames, nodes, dataflowDeps).fixOverlaps().positions();
    setRfNodes(prev => {
      let changed = false;
      const next = prev.map(node => {
        const y = fixed[node.id]?.y;
        if (y === undefined || y === node.position.y) {
          return node;
        }
        changed = true;
        return { ...node, position: { ...node.position, y } };
      });
      return changed ? next : prev;
    });
  }, [measuredKey]);

  // When the canvas height changes, scale the zoom by the same ratio so the same
  // slice of the graph stays framed (anchored at the top-left) instead of
  // revealing more or less of it as the height grows or shrinks.
  useEffect(() => {
    const el = containerRef.current;
    if (!el) {
      return;
    }
    const observer = new ResizeObserver(entries => {
      const height = entries[0].contentRect.height;
      const prev = prevHeightRef.current;
      prevHeightRef.current = height;
      if (prev && height && prev !== height) {
        const ratio = height / prev;
        const { x, y, zoom } = getViewport();
        setViewport({ x: x * ratio, y: y * ratio, zoom: zoom * ratio });
      }
    });
    observer.observe(el);
    return () => observer.disconnect();
  }, [getViewport, setViewport]);

  // Overlay current selection and data (with fresh handlers) each render, so the
  // node always calls the latest rename handler — no stale closures, no ref.
  const displayNodes = rfNodes.map(node => ({
    ...node,
    selected: node.id === selectedNode,
    data: {
      name: node.id,
      node: nodes[node.id] ?? {},
      catalog,
      onRename: (newName: string) => onRenameNode(node.id, newName),
      onChange: (updated: WorkflowNode) => onSetNode(node.id, updated),
      onDelete: () => onDeleteNode(node.id),
    },
  }));

  const rfEdges = useMemo<Edge[]>(() => buildWorkflowEdges(nodeNames, nodes), [structureKey]);

  // Overlay edge type, direction arrow, hover state, and a fresh remove handler.
  const displayEdges = rfEdges.map(edge => ({
    ...edge,
    type: 'requires',
    markerEnd: { type: MarkerType.ArrowClosed },
    data: {
      hovered: edge.id === hoveredEdge,
      onRemove: () => onRemoveRequire(edge.source, edge.target),
    },
  }));

  // Clears the reference a dataflow edge represents from its target parameter —
  // used both by the edge's X button and by deleting the line.
  const removeDataflowEdge = (id: string) => {
    const ref = parseDataflowEdgeId(id);
    if (ref) {
      onSetNode(ref.target, clearInputReference(nodes[ref.target] ?? {}, ref.param, ref.refNode, ref.key));
    }
  };

  // Dataflow edges from parameter values that reference another node's output
  // (e.g. `{C.output.ggg}`), drawn output-handle → input-handle. Computed each
  // render so edits to parameter values re-derive them. Hovering shows an X that
  // clears the reference, mirroring the requires edges.
  const dataflowEdges: Edge[] = dataflowDeps.map(edge => ({
    ...edge,
    type: 'dataflow',
    markerEnd: { type: MarkerType.ArrowClosed, color: '#1976d2' },
    style: { stroke: '#1976d2' },
    animated: true,
    // The input handle sits inside the node, so the line's end runs under the
    // node box; lift it above the nodes so it stays visible.
    zIndex: 1000,
    data: {
      hovered: edge.id === hoveredEdge,
      onRemove: () => removeDataflowEdge(edge.id),
    },
  }));

  const isValidConnection = (connection: Connection | Edge): boolean => {
    // Output→input (dataflow) connections skip the requires cycle check.
    if (parseDataflowConnection(connection.sourceHandle, connection.targetHandle)) {
      return true;
    }
    return isValidConnectionPure(connection, nodeNames, nodes);
  };

  const onConnect = (connection: Connection) => {
    if (!connection.source || !connection.target) {
      return;
    }
    // Dragging an output handle to an input handle writes a dataflow reference
    // ({source.output.name}) into the target's parameter; otherwise it's requires.
    const dataflow = parseDataflowConnection(connection.sourceHandle, connection.targetHandle);
    if (dataflow) {
      onSetNode(connection.target, setInputReference(nodes[connection.target] ?? {}, dataflow.param, connection.source, dataflow.outputName));
      return;
    }
    onAddRequire(connection.source, connection.target);
  };

  const onEdgesDelete = (deleted: Edge[]) => {
    deleted.forEach(edge => {
      // Deleting a dataflow line clears the reference from its parameter; a
      // requires edge removes the requires link.
      if (parseDataflowEdgeId(edge.id)) {
        removeDataflowEdge(edge.id);
        return;
      }
      onRemoveRequire(edge.source, edge.target);
    });
  };

  return (
    <Box ref={containerRef} sx={{ flex: 1, minHeight: 200, ml: -2, mr: -2, mb: -2, borderTop: '1px solid', borderColor: 'divider' }}>
      <ReactFlow
        nodes={displayNodes}
        edges={[...displayEdges, ...dataflowEdges]}
        nodeTypes={NODE_TYPES}
        edgeTypes={EDGE_TYPES}
        onNodesChange={onNodesChange}
        onConnect={onConnect}
        onEdgesDelete={onEdgesDelete}
        isValidConnection={isValidConnection}
        fitView
        onNodeClick={(_e, node) => onSelectNode(node.id)}
        onPaneClick={() => onSelectNode(undefined)}
        onEdgeMouseEnter={(_e, edge) => setHoveredEdge(edge.id)}
        onEdgeMouseLeave={() => setHoveredEdge(null)}
        onNodeContextMenu={(event, node) => {
          event.preventDefault();
          setMenu({ kind: WorkflowContextMenuKind.Node, name: node.id, x: event.clientX, y: event.clientY });
        }}
        onEdgeContextMenu={(event, edge) => {
          event.preventDefault();
          setMenu({ kind: WorkflowContextMenuKind.Edge, source: edge.source, target: edge.target, x: event.clientX, y: event.clientY });
        }}
      >
        <Panel position="top-right">
          <Tooltip title="Add node">
            <IconButton size="small" onClick={onAddNode} sx={{ bgcolor: 'background.paper', boxShadow: 1, mr: 1, p: 0.5 }}>
              <Add fontSize="small" />
            </IconButton>
          </Tooltip>
        </Panel>
        <Background />
        <Controls />
      </ReactFlow>
      <WorkflowContextMenu
        menu={menu}
        onClose={() => setMenu(null)}
        onDeleteNode={onDeleteNode}
        onRemoveRequire={onRemoveRequire}
      />
    </Box>
  );
};

// ReactFlowProvider supplies the store that hooks like useNodesInitialized read,
// so the inner component must live inside it.
export const WorkflowGraph = (props: WorkflowGraphProps) => {
  return (
    <ReactFlowProvider>
      <WorkflowGraphInner {...props} />
    </ReactFlowProvider>
  );
};
