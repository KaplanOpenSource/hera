import { WorkflowNode } from '../../shared/types';
import { NodeCatalogEntry } from './nodeCatalog';
import { Reference } from './references/Reference';
import { ReferenceKind } from './references/ReferenceKind';
import { OUTPUT, knownKinds } from './references/knownKinds';

// A dataflow edge inferred from a parameter value that references another node.
// The reference formats live in the references folder.
export interface WorkflowDataflowEdge {
  id: string;
  source: string;
  sourceHandle: string;
  target: string;
  targetHandle: string;
}

// Handle ids, each scoped to its node so a node's requires and dataflow handles
// never collide. Only the dataflow ids carry a trailing name, so the matchers
// can't mistake a requires handle for a dataflow one.
export const nodeOutputHandleId = (node: string): string => `${node}:req-out`;
export const nodeInputHandleId = (node: string): string => `${node}:req-in`;
export const outputHandleId = (node: string, name: string): string => new Reference(node, OUTPUT, name).handleId();
export const inputHandleId = (node: string, param: string): string => `${node}:in:${param}`;
const INPUT_HANDLE_MATCH = /:in:([^:]+)$/;

// A connection dragged from a source handle to an input handle: the kind and key
// it leaves from and the parameter it lands on.
export interface DataflowConnection {
  kind: ReferenceKind;
  outputName: string;
  param: string;
}

// If a connection runs from a dataflow source handle to a dataflow input handle,
// return the reference/param names it links; otherwise null (it's a requires drag).
export const parseDataflowConnection = (
  sourceHandle: string | null | undefined,
  targetHandle: string | null | undefined,
): DataflowConnection | null => {
  const source = sourceHandle ? knownKinds.ofHandle(sourceHandle) : null;
  const target = targetHandle?.match(INPUT_HANDLE_MATCH);
  if (source && target) {
    return {
      kind: source.kind,
      outputName: source.key,
      param: target[1],
    };
  }
  return null;
};

// The reference token written into an input parameter value to point it at
// another node's output — the same shape buildDataflowEdges parses back out.
// Written as `output`; buildDataflowEdges still accepts the older `parameters`.
export const dataflowReference = (sourceNode: string, outputName: string): string =>
  new Reference(sourceNode, OUTPUT, outputName).toString();

// Splices a reference into `value` at `caret`. The caret is clamped into range,
// so out-of-range positions land at the start or end.
export const insertReferenceAt = (
  value: string,
  caret: number,
  reference: Reference,
): string => {
  const at = Math.max(0, Math.min(caret, value.length));
  return value.slice(0, at) + reference.toString() + value.slice(at);
};

// Which part of a half-typed `{…}` reference the caret sits in: the node name
// (before the first dot), the section, or the key (after the last dot).
export enum ReferenceTokenStage {
  Node = 'node',
  Section = 'section',
  Key = 'key',
}

// The `{…}` reference the caret is inside, as parsed for inline autocomplete.
export interface ReferenceTokenAtCaret {
  stage: ReferenceTokenStage;
  // The node name already typed before the section dot — only set on the Section
  // and Key stages (empty on the Node stage).
  nodePart: string;
  // The section text typed after the node dot — empty on the Node stage.
  sectionPart: string;
  // The kind that section names, or null while no kind reads it yet.
  kind: ReferenceKind | null;
  // The partial text the caret is filtering by: a partial node name (Node stage),
  // a partial section (Section stage) or a partial key (Key stage).
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
    return { stage: ReferenceTokenStage.Node, nodePart: '', sectionPart: '', kind: null, seed: inner, start: open, end };
  }
  const nodePart = inner.slice(0, dot);
  const afterNode = inner.slice(dot + 1);
  // The key is what follows the last dot, so the section may itself hold dots.
  const lastDot = afterNode.lastIndexOf('.');
  const section = lastDot === -1 ? afterNode : afterNode.slice(0, lastDot);
  const kind = knownKinds.bySection(section);
  if (kind === null) {
    return { stage: ReferenceTokenStage.Section, nodePart, sectionPart: afterNode, kind: null, seed: afterNode, start: open, end };
  }
  return {
    stage: ReferenceTokenStage.Key,
    nodePart,
    sectionPart: section,
    kind,
    seed: lastDot === -1 ? '' : afterNode.slice(lastDot + 1),
    start: open,
    end,
  };
};

// Overwrites the token spanning [start, end) with a full reference — used when a
// suggestion is picked from the inline menu.
export const replaceReferenceAt = (
  value: string,
  start: number,
  end: number,
  reference: Reference,
): string => {
  return value.slice(0, start) + reference.toString() + value.slice(end);
};

// Returns node with its `param` input set to the given reference.
export const setInputReference = (
  node: WorkflowNode,
  param: string,
  reference: Reference,
): WorkflowNode => {
  const input_parameters = { ...(node.Execution?.input_parameters ?? {}) };
  input_parameters[param] = reference.toString();
  return { ...node, Execution: { ...node.Execution, input_parameters } };
};

// The parts of a dataflow edge id (df:<refNode>:<mark>:<key>-><target>.<param>).
export interface DataflowEdgeRef {
  refNode: string;
  key: string;
  target: string;
  param: string;
}

// Parses a dataflow edge id back into its parts, or null if it isn't one (e.g. a
// requires edge id) — lets onEdgesDelete tell dataflow lines from requires edges.
export const parseDataflowEdgeId = (id: string): DataflowEdgeRef | null => {
  const parsed = knownKinds.ofEdgeId(id);
  if (parsed === null) {
    return null;
  }
  return { refNode: parsed.reference.node, key: parsed.reference.key, target: parsed.target, param: parsed.param };
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
  const token = new Reference(refNode, OUTPUT, key).clearToken();
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
      for (const reference of knownKinds.parseAll(value)) {
        if (reference.node === target) {
          continue;
        }
        const inGraph = nodeNames.includes(reference.node);
        const isReal = inGraph
          && reference.kind.namesOf(nodes[reference.node] ?? {}, catalog).includes(reference.key);
        const id = reference.edgeIdTo(target, param);
        if (isReal && !seen.has(id)) {
          seen.add(id);
          edges.push({
            id,
            source: reference.node,
            sourceHandle: reference.handleId(),
            target,
            targetHandle: inputHandleId(target, param),
          });
        }
      }
    });
  });
  return edges;
};
