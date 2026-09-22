import { WorkflowNode } from '../../shared/types';
import { insertReferenceAt } from './workflowDataflow';

// Pure edits on a single workflow node. Each one returns a new node; nothing
// here touches the canvas.

// One node's input parameters (never undefined, so callers can spread it).
const paramsOf = (node: WorkflowNode): { [param: string]: any } => {
  return node.Execution?.input_parameters ?? {};
};

// The node with one input parameter dropped (right-click a field -> delete).
export const nodeWithoutParam = (node: WorkflowNode, param: string): WorkflowNode => {
  const input_parameters = { ...paramsOf(node) };
  delete input_parameters[param];
  return { ...node, Execution: { ...node.Execution, input_parameters } };
};

// The node with one input parameter set to a new value.
export const nodeWithParamValue = (node: WorkflowNode, param: string, value: any): WorkflowNode => {
  return { ...node, Execution: { ...node.Execution, input_parameters: { ...paramsOf(node), [param]: value } } };
};

// The node with a {sourceNode.output.name} reference inserted into a parameter
// at the caret (the end when there is none), leaving the rest of the value intact.
export const nodeWithReferenceAt = (
  node: WorkflowNode,
  param: string,
  sourceNode: string,
  output: string,
  caret?: number,
): WorkflowNode => {
  const current = paramsOf(node)[param];
  const text = typeof current === 'string' ? current : '';
  return nodeWithParamValue(node, param, insertReferenceAt(text, caret ?? text.length, sourceNode, output));
};
