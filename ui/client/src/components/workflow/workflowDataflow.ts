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

// Prefix on dataflow edge ids: df:<refNode>.<key>-><target>.<param>.
export const DATAFLOW_EDGE_PREFIX = 'df:';

// Handle ids, each scoped to its node so a node's requires + output/input
// handles never collide. Only the dataflow (:out:/:in:) ids carry a trailing
// name, so the matchers below can't mistake a requires handle for a dataflow one.
export const nodeOutputHandleId = (node: string): string => `${node}:req-out`;
export const nodeInputHandleId = (node: string): string => `${node}:req-in`;
export const outputHandleId = (node: string, name: string): string => `${node}:out:${name}`;
export const inputHandleId = (node: string, param: string): string => `${node}:in:${param}`;
const OUTPUT_HANDLE_MATCH = /:out:([^:]+)$/;
const INPUT_HANDLE_MATCH = /:in:([^:]+)$/;

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
  const source = sourceHandle?.match(OUTPUT_HANDLE_MATCH);
  const target = targetHandle?.match(INPUT_HANDLE_MATCH);
  if (source && target) {
    return {
      outputName: source[1],
      param: target[1],
    };
  }
  return null;
};

// The reference token written into an input parameter value to point it at
// another node's output — the same shape buildDataflowEdges parses back out.
// Written as `parameters` (buildDataflowEdges also accepts `output`).
export const dataflowReference = (sourceNode: string, outputName: string): string =>
  `{${sourceNode}.parameters.${outputName}}`;

// Splices a node-output reference into `value` at `caret`. The caret is clamped
// into range, so out-of-range positions land at the start or end.
export const insertReferenceAt = (
  value: string,
  caret: number,
  sourceNode: string,
  outputName: string,
): string => {
  const token = dataflowReference(sourceNode, outputName);
  const at = Math.max(0, Math.min(caret, value.length));
  return value.slice(0, at) + token + value.slice(at);
};

// Which part of a half-typed `{…}` reference the caret sits in: the node name
// (before the first dot) or the output key (after it).
export enum ReferenceTokenStage {
  Node = 'node',
  Output = 'output',
}

// The `{…}` reference the caret is inside, as parsed for inline autocomplete.
export interface ReferenceTokenAtCaret {
  stage: ReferenceTokenStage;
  // The node name already typed before the section dot — only set on the Output
  // stage (empty on the Node stage).
  nodePart: string;
  // The partial text the caret is filtering by: a partial node name (Node stage)
  // or a partial output key (Output stage).
  seed: string;
  // The token's span in the value, from the opening `{` to just past the closing
  // `}` (or the caret, if the token is still unclosed) — what replaceReferenceAt
  // overwrites when a suggestion is chosen.
  start: number;
  end: number;
}

// If `caret` sits inside a `{…}` reference token, describes it (so the caller can
// show the node or output suggestion menu); otherwise null. A token runs from the
// nearest `{` at/before the caret (with no `}` in between) to the next `}` after
// it, or to the caret itself while still unclosed.
export const tokenAtCaret = (value: string, caret: number): ReferenceTokenAtCaret | null => {
  const at = Math.max(0, Math.min(caret, value.length));
  const open = value.lastIndexOf('{', at - 1);
  if (open === -1 || value.lastIndexOf('}', at - 1) > open) {
    return null;
  }
  const close = value.indexOf('}', open + 1);
  const nextOpen = value.indexOf('{', open + 1);
  const closed = close !== -1 && (nextOpen === -1 || nextOpen > close);
  const end = closed ? close + 1 : at;
  const inner = value.slice(open + 1, at);
  const dot = inner.indexOf('.');
  if (dot === -1) {
    return { stage: ReferenceTokenStage.Node, nodePart: '', seed: inner, start: open, end };
  }
  return {
    stage: ReferenceTokenStage.Output,
    nodePart: inner.slice(0, dot),
    seed: inner.slice(inner.lastIndexOf('.') + 1),
    start: open,
    end,
  };
};

// Overwrites the token spanning [start, end) with a full reference to sourceNode's
// output — used when a suggestion is picked from the inline menu.
export const replaceReferenceAt = (
  value: string,
  start: number,
  end: number,
  sourceNode: string,
  outputName: string,
): string => {
  return value.slice(0, start) + dataflowReference(sourceNode, outputName) + value.slice(end);
};

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
            sourceHandle: outputHandleId(refNode, key),
            target,
            targetHandle: inputHandleId(target, param),
          });
        }
        match = REFERENCE.exec(value);
      }
    });
  });
  return edges;
};
