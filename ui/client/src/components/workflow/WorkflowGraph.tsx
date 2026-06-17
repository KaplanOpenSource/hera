import { Box } from '@mui/material';
import { Background, Controls, Edge, Node, ReactFlow, useNodesState } from '@xyflow/react';
import '@xyflow/react/dist/style.css';
import { useEffect, useMemo } from 'react';
import { WorkflowNode } from '../../shared/types';

const X_GAP = 240;
const Y_GAP = 90;

const normalizeRequires = (requires?: string | string[]): string[] => {
  if (requires === undefined) {
    return [];
  }
  return Array.isArray(requires) ? requires : [requires];
};

// Assigns each node a layer = longest `requires` chain depth, so the graph lays
// out left-to-right by dependency order. Cycles are broken at layer 0.
const computeLayers = (nodeNames: string[], nodes: Record<string, WorkflowNode>): Record<string, number> => {
  const layer: Record<string, number> = {};
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
// `requires` field. Nodes are draggable; clicking a node selects it for editing.
export const WorkflowGraph = ({
  nodeNames,
  nodes,
  selectedNode,
  onSelectNode,
}: {
  nodeNames: string[],
  nodes: Record<string, WorkflowNode>,
  selectedNode?: string,
  onSelectNode: (name: string) => void,
}) => {
  const [rfNodes, setRfNodes, onNodesChange] = useNodesState<Node>([]);

  // A signature of the workflow structure (names, types, requires) so the graph
  // only rebuilds when the structure changes — not when a node is dragged.
  const structureKey = JSON.stringify(nodeNames.map(name => [name, nodes[name]?.type, nodes[name]?.requires]));

  // Rebuild nodes on structure change, keeping positions the user has dragged.
  useEffect(() => {
    const layers = computeLayers(nodeNames, nodes);
    const rowInLayer: Record<number, number> = {};
    setRfNodes(prev => {
      const prevPos = new Map(prev.map(n => [n.id, n.position]));
      return nodeNames.map(name => {
        const layer = layers[name] ?? 0;
        const row = rowInLayer[layer] ?? 0;
        rowInLayer[layer] = row + 1;
        return {
          id: name,
          position: prevPos.get(name) ?? { x: layer * X_GAP, y: row * Y_GAP },
          selected: name === selectedNode,
          data: {
            label: (
              <div>
                {name}
                <div style={{ fontSize: 10, opacity: 0.7 }}>{nodes[name]?.type || '—'}</div>
              </div>
            ),
          },
        };
      });
    });
  }, [structureKey]);

  // Reflect selection without rebuilding (so it doesn't reset dragged positions).
  useEffect(() => {
    setRfNodes(prev => prev.map(n => ({ ...n, selected: n.id === selectedNode })));
  }, [selectedNode]);

  const rfEdges = useMemo<Edge[]>(() => {
    const edges: Edge[] = [];
    nodeNames.forEach(name => {
      normalizeRequires(nodes[name]?.requires)
        .filter(req => nodeNames.includes(req))
        .forEach(req => edges.push({ id: `${req}->${name}`, source: req, target: name }));
    });
    return edges;
  }, [structureKey]);

  return (
    <Box sx={{ height: 400, border: '1px solid', borderColor: 'divider', mb: 2 }}>
      <ReactFlow
        nodes={rfNodes}
        edges={rfEdges}
        onNodesChange={onNodesChange}
        fitView
        onNodeClick={(_e, node) => onSelectNode(node.id)}
      >
        <Background />
        <Controls />
      </ReactFlow>
    </Box>
  );
};
