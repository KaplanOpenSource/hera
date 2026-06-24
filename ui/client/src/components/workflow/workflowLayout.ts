import { WorkflowNode } from '../../shared/types';
import { normalizeRequires } from '../../shared/workflow';

export const X_GAP = 700;       // horizontal distance between dependency layers
export const V_GAP = 100;       // vertical gap between nodes in a column
export const BASE_HEIGHT = 110; // node height without params (name + type)
export const ROW_HEIGHT = 28;   // estimated height per parameter row

// Assigns each node a layer = longest `requires` chain depth, so the graph lays
// out left-to-right by dependency order. Cycles are broken at layer 0.
export const computeLayers = (nodeNames: string[], nodes: { [name: string]: WorkflowNode }): { [name: string]: number } => {
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

// Number of rows a parameter tree renders to (one per key, recursively), used to
// estimate a node's height before it is measured (avoids a layout flash).
export const countRows = (value: unknown): number => {
  if (value && typeof value === 'object') {
    return Object.values(value).reduce((sum: number, child) => sum + 1 + countRows(child), 0);
  }
  return 0;
};

// Estimated render height of a node from its parameter tree.
export const estimateHeight = (node: WorkflowNode): number => {
  const params = node.Execution?.input_parameters ?? {};
  return BASE_HEIGHT + (1 + countRows(params)) * ROW_HEIGHT;
};

// Stacks each dependency layer into a column using estimated node heights.
export const computeLayout = (nodeNames: string[], nodes: { [name: string]: WorkflowNode }): { [name: string]: { x: number, y: number } } => {
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

