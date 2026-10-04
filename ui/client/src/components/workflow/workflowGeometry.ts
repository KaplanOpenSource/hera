import { WorkflowNode } from '../../shared/types';
import { normalizeRequires } from '../../shared/workflow';

// Sized for the summary card a node shows when the pointer is away from it, not
// for the editor that opens over it on hover - the editor is laid over the
// canvas and needs no room of its own.
export const X_GAP = 420;      // horizontal distance between dependency layers
export const V_GAP = 40;       // vertical gap between nodes in a column
export const BASE_HEIGHT = 52; // card height without params (name + type)
export const ROW_HEIGHT = 22;  // card height per parameter row

// Assigns each node a layer = longest dependency chain depth, so the graph lays
// out left-to-right by dependency order. Dependencies are `requires` plus any
// extraDeps (e.g. dataflow references). Cycles are broken at layer 0.
export const computeLayers = (
  nodeNames: string[],
  nodes: { [name: string]: WorkflowNode },
  extraDeps: { source: string, target: string }[] = [],
): { [name: string]: number } => {
  const extraPreds: { [target: string]: string[] } = {};
  extraDeps.forEach(({ source, target }) => {
    (extraPreds[target] ??= []).push(source);
  });
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
    const reqs = [...normalizeRequires(nodes[name]?.requires), ...(extraPreds[name] ?? [])]
      .filter(r => nodeNames.includes(r) && r !== name);
    const value = reqs.length === 0 ? 0 : Math.max(...reqs.map(resolve)) + 1;
    visiting.delete(name);
    layer[name] = value;
    return value;
  };
  nodeNames.forEach(resolve);
  return layer;
};

// Estimated height of a node's summary card before it is measured (avoids a
// layout flash). The card lists one row per parameter, nested values and all
// their rows staying hidden until the editor opens.
export const estimateHeight = (node: WorkflowNode): number => {
  const params = node.Execution?.input_parameters ?? {};
  return BASE_HEIGHT + Object.keys(params).length * ROW_HEIGHT;
};

// WorkflowLayout (in ./WorkflowLayout) builds on these helpers to place nodes and
// resolve overlaps.

