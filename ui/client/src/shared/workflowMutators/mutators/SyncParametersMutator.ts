import { WorkflowBlock, WorkflowDesc } from '../../types';
import { getWorkflowBlock } from '../../workflow';

// Each node's input_parameters keyed by node name. Mirrors the Hermes
// parametersJSON index the database stores under desc.parameters.
export const workflowParameters = (block: WorkflowBlock): { [node: string]: any } => {
  const nodes = block.nodes ?? {};
  const names = block.nodeList ?? Object.keys(nodes);
  const parameters: { [node: string]: any } = {};
  for (const name of names) {
    const params = nodes[name]?.Execution?.input_parameters;
    if (params !== undefined) {
      parameters[name] = params;
    }
  }
  return parameters;
};

// Rebuilds desc.parameters from the workflow's nodes, so it never drifts. A
// desc with no block is left alone.
export const syncParameters = (desc: WorkflowDesc): WorkflowDesc => {
  const block = getWorkflowBlock(desc.workflow);
  if (!block) {
    return desc;
  }
  return { ...desc, parameters: workflowParameters(block) };
};
