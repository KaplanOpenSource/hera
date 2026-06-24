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
import { WorkflowLayout } from './WorkflowLayout';

// Defined once (module scope) so ReactFlow doesn't warn about changing types.
const NODE_TYPES = { workflow: WorkflowFlowNode };
const EDGE_TYPES = { requires: WorkflowRequiresEdge };

interface WorkflowGraphProps {
  catalog: NodeCatalogEntry[];
  nodeNames: string[];
  nodes: { [name: string]: WorkflowNode };
  selectedNode?: string;
  onSelectNode: (name: string) => void;
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

  // A signature of the workflow structure (names, types, requires) so the graph
  // only rebuilds when the structure changes — not when a node is dragged.
  const structureKey = JSON.stringify(nodeNames.map(name => [name, nodes[name]?.type, nodes[name]?.requires]));

  useEffect(() => {
    const layout = WorkflowLayout.stacked(nodeNames, nodes).positions();
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
  const membershipKey = JSON.stringify([...nodeNames].sort());
  useEffect(() => {
    const layout = WorkflowLayout.stacked(nodeNames, nodes).positions();
    setRfNodes(prev => prev.map(node => ({ ...node, position: layout[node.id] ?? node.position })));
    const isInitial = prevNamesRef.current.length === 0;
    const added = nodeNames.filter(n => !prevNamesRef.current.includes(n));
    const removed = prevNamesRef.current.filter(n => !nodeNames.includes(n));
    prevNamesRef.current = nodeNames;
    if (isInitial) {
      pendingRef.current = 'all';
    } else if (added.length > 0 && removed.length === 0) {
      pendingRef.current = added.length === 1 ? added[0] : 'all';
    } else if (removed.length > 0 && added.length === 0) {
      requestAnimationFrame(() => fitView({ duration: 300 }));
    }
  }, [membershipKey]);

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
    const fixed = WorkflowLayout.fromFlowNodes(rfNodes, nodeNames, nodes).fixOverlaps().positions();
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

  const isValidConnection = (connection: Connection | Edge): boolean => {
    return isValidConnectionPure(connection, nodeNames, nodes);
  };

  const onConnect = (connection: Connection) => {
    if (connection.source && connection.target) {
      onAddRequire(connection.source, connection.target);
    }
  };

  const onEdgesDelete = (deleted: Edge[]) => {
    deleted.forEach(edge => onRemoveRequire(edge.source, edge.target));
  };

  return (
    <Box ref={containerRef} sx={{ flex: 1, minHeight: 200, ml: -2, mr: -2, mb: -2, borderTop: '1px solid', borderColor: 'divider' }}>
      <ReactFlow
        nodes={displayNodes}
        edges={displayEdges}
        nodeTypes={NODE_TYPES}
        edgeTypes={EDGE_TYPES}
        onNodesChange={onNodesChange}
        onConnect={onConnect}
        onEdgesDelete={onEdgesDelete}
        isValidConnection={isValidConnection}
        fitView
        onNodeClick={(_e, node) => onSelectNode(node.id)}
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
