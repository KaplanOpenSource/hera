import { Add } from '@mui/icons-material';
import { Box, IconButton, Tooltip } from '@mui/material';
import { Background, Connection, Controls, Edge, Node, Panel, ReactFlow, useNodesState } from '@xyflow/react';
import '@xyflow/react/dist/style.css';
import { useEffect, useMemo } from 'react';
import { WorkflowNode } from '../../shared/types';
import { normalizeRequires } from '../../shared/workflow';
import { WorkflowFlowNode } from './WorkflowFlowNode';

const X_GAP = 240;
const Y_GAP = 90;

// Defined once (module scope) so ReactFlow doesn't warn about changing nodeTypes.
const NODE_TYPES = { workflow: WorkflowFlowNode };

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

// Node graph view of a workflow: one node per workflow node, edges from the
// `requires` field. Nodes are draggable, their names editable inline, and the
// on-canvas button adds a node. Clicking a node selects it for editing.
export const WorkflowGraph = ({
  nodeNames,
  nodes,
  selectedNode,
  onSelectNode,
  onAddNode,
  onRenameNode,
  onAddRequire,
  onRemoveRequire,
}: {
  nodeNames: string[],
  nodes: { [name: string]: WorkflowNode },
  selectedNode?: string,
  onSelectNode: (name: string) => void,
  onAddNode: () => void,
  onRenameNode: (oldName: string, newName: string) => void,
  onAddRequire: (source: string, target: string) => void,
  onRemoveRequire: (source: string, target: string) => void,
}) => {
  // useNodesState owns only position and identity; the structure effect rebuilds
  // it (preserving dragged positions) when nodes are added/removed/reordered.
  const [rfNodes, setRfNodes, onNodesChange] = useNodesState<Node>([]);

  // A signature of the workflow structure (names, types, requires) so the graph
  // only rebuilds when the structure changes — not when a node is dragged.
  const structureKey = JSON.stringify(nodeNames.map(name => [name, nodes[name]?.type, nodes[name]?.requires]));

  useEffect(() => {
    const layers = computeLayers(nodeNames, nodes);
    const rowInLayer: { [layer: number]: number } = {};
    setRfNodes(prev => {
      const prevPos = new Map(prev.map(n => [n.id, n.position]));
      return nodeNames.map(name => {
        const layer = layers[name] ?? 0;
        const row = rowInLayer[layer] ?? 0;
        rowInLayer[layer] = row + 1;
        return {
          id: name,
          type: 'workflow',
          position: prevPos.get(name) ?? { x: layer * X_GAP, y: row * Y_GAP },
          data: {},
        };
      });
    });
  }, [structureKey]);

  // Overlay current selection and data (with fresh handlers) each render, so the
  // node always calls the latest rename handler — no stale closures, no ref.
  const displayNodes = rfNodes.map(node => ({
    ...node,
    selected: node.id === selectedNode,
    data: {
      name: node.id,
      type: nodes[node.id]?.type,
      onRename: (newName: string) => onRenameNode(node.id, newName),
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
        edges={rfEdges}
        nodeTypes={NODE_TYPES}
        onNodesChange={onNodesChange}
        onConnect={onConnect}
        onEdgesDelete={onEdgesDelete}
        isValidConnection={isValidConnection}
        fitView
        onNodeClick={(_e, node) => onSelectNode(node.id)}
      >
        <Panel position="top-right">
          <Tooltip title="Add node">
            <IconButton size="small" onClick={onAddNode} sx={{ bgcolor: 'background.paper', boxShadow: 1 }}>
              <Add />
            </IconButton>
          </Tooltip>
        </Panel>
        <Background />
        <Controls />
      </ReactFlow>
    </Box>
  );
};
