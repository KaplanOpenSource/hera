import { WorkflowNode } from '../../shared/types';
import { normalizeRequires } from '../../shared/workflow';

// A graph edge derived from a node's requires (source must run before target).
export interface WorkflowEdge {
  id: string;
  source: string;
  target: string;
}

// Builds one edge per requires link (source→target = "target requires source"),
// skipping requires that point to nodes outside the graph.
export const buildWorkflowEdges = (nodeNames: string[], nodes: { [name: string]: WorkflowNode }): WorkflowEdge[] => {
  const edges: WorkflowEdge[] = [];
  nodeNames.forEach(name => {
    normalizeRequires(nodes[name]?.requires)
      .filter(req => nodeNames.includes(req))
      .forEach(req => edges.push({ id: `${req}->${name}`, source: req, target: name }));
  });
  return edges;
};

// A connection source→target means "target requires source". Reject it if it
// points the wrong way — i.e. would create a cycle because target can already
// reach source through existing requires — or is a self/duplicate edge.
export const isValidConnection = (
  connection: { source?: string | null, target?: string | null },
  nodeNames: string[],
  nodes: { [name: string]: WorkflowNode },
): boolean => {
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
