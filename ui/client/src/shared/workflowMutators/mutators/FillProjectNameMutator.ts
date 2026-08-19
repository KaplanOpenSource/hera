import { WorkflowBlock, WorkflowDesc } from '../../types';
import { WorkflowMutatorBase, WorkflowMutatorContext } from '../WorkflowMutatorBase';

// True for a parameter key that names a project: `ProjectName`, `projectName`,
// `projectname` - any casing. Hera/Hermes nodes use `ProjectName`.
export const isProjectNameKey = (key: string): boolean => {
  return /^projectname$/i.test(key);
};

// Sets every project-name parameter (any casing) in the block to the given
// project. Looks at each node's top-level input_parameters and overwrites any
// existing value.
export const fillProjectName = (block: WorkflowBlock, projectName: string): WorkflowBlock => {
  const nodes = { ...block.nodes };
  for (const [name, node] of Object.entries(nodes)) {
    const params = node.Execution?.input_parameters;
    if (!params) {
      continue;
    }
    const nextParams = { ...params };
    let changed = false;
    for (const key of Object.keys(nextParams)) {
      if (isProjectNameKey(key)) {
        nextParams[key] = projectName;
        changed = true;
      }
    }
    if (changed) {
      nodes[name] = {
        ...node,
        Execution: { ...node.Execution, input_parameters: nextParams },
      };
    }
  }
  return { ...block, nodes };
};

// Phase: set every project-name param (any casing) to the project the document
// lives in, so the field always matches its project.
export class FillProjectNameMutator extends WorkflowMutatorBase {
  readonly name = 'fillProjectName';

  mutate(desc: WorkflowDesc, context: WorkflowMutatorContext): WorkflowDesc {
    const block = this.getBlock(desc);
    if (!block || !context.projectName) {
      return desc;
    }
    return this.writeBlock(desc, fillProjectName(block, context.projectName));
  }
}
