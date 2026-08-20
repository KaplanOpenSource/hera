import { WorkflowBlock, WorkflowDesc } from '../../types';
import { WorkflowMutatorBase, WorkflowMutatorContext } from '../WorkflowMutatorBase';

// The `parameters` map derived from a block: each node's input_parameters keyed
// by node name. Mirrors the Hermes `workflow.parametersJSON` property, which the
// database stores under `desc.parameters` and queries workflows by. Kept in sync
// so a query on parameters matches the workflow's real, current node params.
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

// Phase: rebuild desc.parameters (the query index) from the workflow's nodes, so
// it never drifts from the workflow itself.
export class SyncParametersMutator extends WorkflowMutatorBase {
  readonly name = 'syncParameters';

  // Ignores context: the parameters index depends only on the workflow itself.
  mutate(desc: WorkflowDesc, _context: WorkflowMutatorContext): WorkflowDesc {
    const block = this.getBlock(desc);
    if (!block) {
      return desc;
    }
    return { ...desc, parameters: workflowParameters(block) };
  }
}
