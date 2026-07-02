import { WorkflowNode } from '../../shared/types';
import { NodeCatalogEntry, nodeOutputNames } from './nodeCatalog';

// A dataflow edge inferred from an input parameter that references another
// node's output — e.g. a param value containing `{C.output.ggg}` (also
// `parameter(s)` / `outputs`), where node C produces an output named `ggg`.
export interface WorkflowDataflowEdge {
  id: string;
  source: string;
  sourceHandle: string;
  target: string;
  targetHandle: string;
}

// Handle ids (must match the ids rendered on the output chips / input rows).
export const outputHandleId = (name: string): string => `out:${name}`;
export const inputHandleId = (name: string): string => `in:${name}`;

// A reference embedded in a parameter value: {<node>.<section>.<key>}, where the
// section is parameter(s) or output(s). We only care about the node and the key.
const REFERENCE = /\{\s*(\w+)\.(?:parameters?|outputs?)\.(\w+)\s*\}/g;

// Builds one edge per input parameter whose value references an output of another
// node in the graph. Only top-level string parameter values are scanned.
export const buildDataflowEdges = (
  nodeNames: string[],
  nodes: { [name: string]: WorkflowNode },
  catalog: NodeCatalogEntry[],
): WorkflowDataflowEdge[] => {
  const edges: WorkflowDataflowEdge[] = [];
  const seen = new Set<string>();
  nodeNames.forEach(target => {
    const params = nodes[target]?.Execution?.input_parameters ?? {};
    Object.entries(params).forEach(([param, value]) => {
      if (typeof value !== 'string') {
        return;
      }
      REFERENCE.lastIndex = 0;
      let match = REFERENCE.exec(value);
      while (match !== null) {
        const refNode = match[1];
        const key = match[2];
        const hasOutput = nodeNames.includes(refNode) && nodeOutputNames(nodes[refNode] ?? {}, catalog).includes(key);
        const id = `df:${refNode}.${key}->${target}.${param}`;
        if (hasOutput && !seen.has(id)) {
          seen.add(id);
          edges.push({
            id,
            source: refNode,
            sourceHandle: outputHandleId(key),
            target,
            targetHandle: inputHandleId(param),
          });
        }
        match = REFERENCE.exec(value);
      }
    });
  });
  return edges;
};
