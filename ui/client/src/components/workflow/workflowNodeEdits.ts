import { WorkflowNode } from '../../shared/types';
import { insertReferenceAt } from './workflowDataflow';
import { Reference } from './references/Reference';
import { knownKinds } from './references/knownKinds';
import { valueAtPath, withoutPath, withValueAtPath } from './paramPath';
import { normalizeRequires } from '../../shared/workflow';

// Pure edits on a single workflow node. Each one returns a new node; nothing
// here touches the canvas. A `paramPath` names a top-level parameter, or a
// value nested inside one (see paramPath.ts).

// One node's input parameters (never undefined, so callers can spread it).
const paramsOf = (node: WorkflowNode): { [param: string]: any } => {
  return node.Execution?.input_parameters ?? {};
};

// The node with one input value dropped (right-click a field -> delete).
export const nodeWithoutParam = (node: WorkflowNode, paramPath: string): WorkflowNode => {
  return { ...node, Execution: { ...node.Execution, input_parameters: withoutPath(paramsOf(node), paramPath) } };
};

// The node with one input value set to a new value.
export const nodeWithParamValue = (node: WorkflowNode, paramPath: string, value: any): WorkflowNode => {
  return { ...node, Execution: { ...node.Execution, input_parameters: withValueAtPath(paramsOf(node), paramPath, value) } };
};

// The node with a reference inserted into a value at the caret (the end when
// there is none), leaving the rest of the value intact.
export const nodeWithReferenceAt = (
  node: WorkflowNode,
  paramPath: string,
  reference: Reference,
  caret?: number,
): WorkflowNode => {
  const current = valueAtPath(paramsOf(node), paramPath);
  const text = typeof current === 'string' ? current : '';
  return nodeWithParamValue(node, paramPath, insertReferenceAt(text, caret ?? text.length, reference));
};

// The node with one more name in `requires`, unchanged when it is already there
// or is the node itself. A reference to another node's output also means the
// other node must run first, so writing a reference goes through here too.
export const nodeWithRequire = (node: WorkflowNode, name: string, selfName: string): WorkflowNode => {
  const current = normalizeRequires(node.requires);
  if (name === selfName || current.includes(name)) {
    return node;
  }
  return { ...node, requires: [...current, name] };
};

// The node with a requires link for every node a text value references, so a
// typed or picked `{other.output}` means "other runs first" as well. Only names
// that are nodes of this workflow count.
export const nodeWithRequiresFromText = (
  node: WorkflowNode,
  selfName: string,
  text: string,
  nodeNames: string[],
): WorkflowNode => {
  let next = node;
  for (const reference of knownKinds.parseAll(text)) {
    if (nodeNames.includes(reference.node)) {
      next = nodeWithRequire(next, reference.node, selfName);
    }
  }
  return next;
};
