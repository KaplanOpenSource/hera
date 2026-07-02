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

// Prefixes marking a dataflow handle id (as opposed to a node-level requires
// handle, which has no id). Output handles leave a node; input handles land.
export const OUTPUT_HANDLE_PREFIX = 'out:';
export const INPUT_HANDLE_PREFIX = 'in:';

// Prefix on dataflow edge ids: df:<refNode>.<key>-><target>.<param>.
export const DATAFLOW_EDGE_PREFIX = 'df:';

// Handle ids (must match the ids rendered on the output chips / input rows).
export const outputHandleId = (name: string): string => `${OUTPUT_HANDLE_PREFIX}${name}`;
export const inputHandleId = (name: string): string => `${INPUT_HANDLE_PREFIX}${name}`;

// A reference embedded in a parameter value: {<node>.<section>.<key>}, where the
// section is parameter(s) or output(s). We only care about the node and the key.
const REFERENCE = /\{\s*(\w+)\.(?:parameters?|outputs?)\.(\w+)\s*\}/g;

// A connection dragged from an output handle to an input handle: the output name
// it leaves from and the parameter it lands on.
export interface DataflowConnection {
  outputName: string;
  param: string;
}

// If a connection runs from a dataflow output handle to a dataflow input handle,
// return the output/param names it links; otherwise null (it's a requires drag).
export const parseDataflowConnection = (
  sourceHandle: string | null | undefined,
  targetHandle: string | null | undefined,
): DataflowConnection | null => {
  if (sourceHandle?.startsWith(OUTPUT_HANDLE_PREFIX) && targetHandle?.startsWith(INPUT_HANDLE_PREFIX)) {
    return {
      outputName: sourceHandle.slice(OUTPUT_HANDLE_PREFIX.length),
      param: targetHandle.slice(INPUT_HANDLE_PREFIX.length),
    };
  }
  return null;
};

// The reference token written into an input parameter value to point it at
// another node's output — the same shape buildDataflowEdges parses back out.
// Written as `parameters` (buildDataflowEdges also accepts `output`).
export const dataflowReference = (sourceNode: string, outputName: string): string =>
  `{${sourceNode}.parameters.${outputName}}`;

// Returns node with its `param` input set to reference sourceNode's output.
export const setInputReference = (
  node: WorkflowNode,
  param: string,
  sourceNode: string,
  outputName: string,
): WorkflowNode => {
  const input_parameters = { ...(node.Execution?.input_parameters ?? {}) };
  input_parameters[param] = dataflowReference(sourceNode, outputName);
  return { ...node, Execution: { ...node.Execution, input_parameters } };
};

// The parts of a dataflow edge id (df:<refNode>.<key>-><target>.<param>).
export interface DataflowEdgeRef {
  refNode: string;
  key: string;
  target: string;
  param: string;
}

// Parses a dataflow edge id back into its parts, or null if it isn't one (e.g. a
// requires edge id) — lets onEdgesDelete tell dataflow lines from requires edges.
export const parseDataflowEdgeId = (id: string): DataflowEdgeRef | null => {
  if (!id.startsWith(DATAFLOW_EDGE_PREFIX)) {
    return null;
  }
  const match = /^df:(\w+)\.(\w+)->(\w+)\.(\w+)$/.exec(id);
  if (match === null) {
    return null;
  }
  return { refNode: match[1], key: match[2], target: match[3], param: match[4] };
};

// Returns node with the reference to refNode's `key` removed from `param`'s value
// — the inverse of setInputReference when a dataflow line is deleted. Only the
// matching {refNode.(parameters|outputs).key} token is stripped from the string.
export const clearInputReference = (
  node: WorkflowNode,
  param: string,
  refNode: string,
  key: string,
): WorkflowNode => {
  const params = node.Execution?.input_parameters ?? {};
  const value = params[param];
  if (typeof value !== 'string') {
    return node;
  }
  const token = new RegExp(`\\{\\s*${refNode}\\.(?:parameters?|outputs?)\\.${key}\\s*\\}`, 'g');
  const input_parameters = { ...params, [param]: value.replace(token, '').trim() };
  return { ...node, Execution: { ...node.Execution, input_parameters } };
};

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
        const id = `${DATAFLOW_EDGE_PREFIX}${refNode}.${key}->${target}.${param}`;
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
