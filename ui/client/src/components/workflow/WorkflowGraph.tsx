import { Add, AutoFixHigh } from '@mui/icons-material';
import { Box, IconButton, Menu, MenuItem, Tooltip } from '@mui/material';
import { Background, Connection, Controls, Edge, MarkerType, Node, Panel, ReactFlow, ReactFlowProvider, useNodesState, useReactFlow } from '@xyflow/react';
import '@xyflow/react/dist/style.css';
import { useEffect, useMemo, useState } from 'react';
import { WorkflowNode } from '../../shared/types';
import { normalizeRequires } from '../../shared/workflow';
import { WorkflowFlowNode } from './WorkflowFlowNode';
import { WorkflowRequiresEdge } from './WorkflowRequiresEdge';

const X_GAP = 340;       // horizontal distance between dependency layers
const V_GAP = 30;        // vertical gap between nodes in a column
const BASE_HEIGHT = 110; // node height without params (name + type)
const ROW_HEIGHT = 28;   // estimated height per parameter row

// Defined once (module scope) so ReactFlow doesn't warn about changing types.
const NODE_TYPES = { workflow: WorkflowFlowNode };
const EDGE_TYPES = { requires: WorkflowRequiresEdge };

// Right-click target: a node, or an edge (a requires link), anchored at a point.
type ContextMenu =
  | { kind: 'node', name: string, x: number, y: number }
  | { kind: 'edge', source: string, target: string, x: number, y: number };

// Assigns each node a layer = longest `requires` chain depth, so the graph lays
// out left-to-right by dependency order. Cycles are broken at layer 0.
const computeLayers = (nodeNames: string[], nodes: { [name: string]: WorkflowNode }): { [name: string]: number } => {
  const layer: { [name: string]: number } = {};
  const visiting = new Set<string>();
  const resolve = (name: string): number => {
    if (layer[name] !== undefined) {
      return layer[name];
    }
    if (visiting.has(name)) {
      return 0;
    }
    visiting.add(name);
    const reqs = normalizeRequires(nodes[name]?.requires).filter(r => nodeNames.includes(r));
    const value = reqs.length === 0 ? 0 : Math.max(...reqs.map(resolve)) + 1;
    visiting.delete(name);
    layer[name] = value;
    return value;
  };
  nodeNames.forEach(resolve);
  return layer;
};

// Estimated render height of a node from its parameter tree, so columns can be
// laid out before the nodes are measured (avoids a layout flash).
const countRows = (value: unknown): number => {
  if (value && typeof value === 'object') {
    return Object.values(value).reduce((sum: number, child) => sum + 1 + countRows(child), 0);
  }
  return 0;
};

const estimateHeight = (node: WorkflowNode): number => {
  const params = node.Execution?.input_parameters ?? {};
  return BASE_HEIGHT + (1 + countRows(params)) * ROW_HEIGHT;
};

// Stacks each dependency layer into a column using estimated node heights.
const computeLayout = (nodeNames: string[], nodes: { [name: string]: WorkflowNode }): { [name: string]: { x: number, y: number } } => {
  const layers = computeLayers(nodeNames, nodes);
  const columnBottom: { [layer: number]: number } = {};
  const positions: { [name: string]: { x: number, y: number } } = {};
  nodeNames.forEach(name => {
    const layer = layers[name] ?? 0;
    const y = columnBottom[layer] ?? 0;
    positions[name] = { x: layer * X_GAP, y };
    columnBottom[layer] = y + estimateHeight(nodes[name] ?? {}) + V_GAP;
  });
  return positions;
};

interface WorkflowGraphProps {
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
  const [menu, setMenu] = useState<ContextMenu | null>(null);
  const [hoveredEdge, setHoveredEdge] = useState<string | null>(null);
  const { fitView } = useReactFlow();
  // useNodesState owns only position and identity; the structure effect rebuilds
  // it (preserving dragged positions) when nodes are added/removed/reordered.
  const [rfNodes, setRfNodes, onNodesChange] = useNodesState<Node>([]);

  // A signature of the workflow structure (names, types, requires) so the graph
  // only rebuilds when the structure changes — not when a node is dragged.
  const structureKey = JSON.stringify(nodeNames.map(name => [name, nodes[name]?.type, nodes[name]?.requires]));

  useEffect(() => {
    const layout = computeLayout(nodeNames, nodes);
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

  // Re-stack all nodes by estimated height and zoom to fit.
  const tidyLayout = () => {
    const layout = computeLayout(nodeNames, nodes);
    setRfNodes(prev => prev.map(node => ({ ...node, position: layout[node.id] ?? node.position })));
    requestAnimationFrame(() => fitView({ duration: 300 }));
  };

  // Tidy and zoom to fit whenever nodes are added or removed (and on mount) —
  // but not on every drag or param edit.
  const membershipKey = JSON.stringify([...nodeNames].sort());
  useEffect(() => {
    tidyLayout();
  }, [membershipKey]);

  // Overlay current selection and data (with fresh handlers) each render, so the
  // node always calls the latest rename handler — no stale closures, no ref.
  const displayNodes = rfNodes.map(node => ({
    ...node,
    selected: node.id === selectedNode,
    data: {
      name: node.id,
      node: nodes[node.id] ?? {},
      onRename: (newName: string) => onRenameNode(node.id, newName),
      onChange: (updated: WorkflowNode) => onSetNode(node.id, updated),
      onDelete: () => onDeleteNode(node.id),
    },
  }));

  const rfEdges = useMemo<Edge[]>(() => {
    const edges: Edge[] = [];
    nodeNames.forEach(name => {
      normalizeRequires(nodes[name]?.requires)
        .filter(req => nodeNames.includes(req))
        .forEach(req => edges.push({ id: `${req}->${name}`, source: req, target: name }));
    });
    return edges;
  }, [structureKey]);

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

  // A connection source→target means "target requires source". Reject it if it
  // points the wrong way — i.e. would create a cycle because target can already
  // reach source through existing requires.
  const isValidConnection = (connection: Connection | Edge): boolean => {
    const { source, target } = connection;
    if (!source || !target || source === target) {
      return false;
    }
    if (normalizeRequires(nodes[target]?.requires).includes(source)) {
      return false;
    }
    const successors: { [name: string]: string[] } = {};
    nodeNames.forEach(name => {
      normalizeRequires(nodes[name]?.requires).forEach(pred => {
        if (!successors[pred]) {
          successors[pred] = [];
        }
        successors[pred].push(name);
      });
    });
    const reaches = (from: string, goal: string): boolean => {
      const seen = new Set<string>();
      const stack = [from];
      while (stack.length > 0) {
        const current = stack.pop() as string;
        if (current === goal) {
          return true;
        }
        if (seen.has(current)) {
          continue;
        }
        seen.add(current);
        (successors[current] ?? []).forEach(next => stack.push(next));
      }
      return false;
    };
    return !reaches(target, source);
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
    <Box sx={{ height: 400, border: '1px solid', borderColor: 'divider', mb: 2 }}>
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
          setMenu({ kind: 'node', name: node.id, x: event.clientX, y: event.clientY });
        }}
        onEdgeContextMenu={(event, edge) => {
          event.preventDefault();
          setMenu({ kind: 'edge', source: edge.source, target: edge.target, x: event.clientX, y: event.clientY });
        }}
      >
        <Panel position="top-right">
          <Tooltip title="Add node">
            <IconButton size="small" onClick={onAddNode} sx={{ bgcolor: 'background.paper', boxShadow: 1, mr: 1 }}>
              <Add />
            </IconButton>
          </Tooltip>
          {/* <Tooltip title="Tidy layout">
            <IconButton size="small" onClick={tidyLayout} sx={{ bgcolor: 'background.paper', boxShadow: 1 }}>
              <AutoFixHigh />
            </IconButton>
          </Tooltip> */}
        </Panel>
        <Background />
        <Controls />
      </ReactFlow>
      <Menu
        open={menu !== null}
        onClose={() => setMenu(null)}
        anchorReference="anchorPosition"
        anchorPosition={menu ? { top: menu.y, left: menu.x } : undefined}
      >
        {menu?.kind === 'node' && (
          <MenuItem onClick={() => { onDeleteNode(menu.name); setMenu(null); }}>
            Delete node “{menu.name}”
          </MenuItem>
        )}
        {menu?.kind === 'edge' && (
          <MenuItem onClick={() => { onRemoveRequire(menu.source, menu.target); setMenu(null); }}>
            Remove requirement ({menu.source} → {menu.target})
          </MenuItem>
        )}
      </Menu>
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
