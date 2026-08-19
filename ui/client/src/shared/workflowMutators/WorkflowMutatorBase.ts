import { WorkflowBlock, WorkflowDesc } from '../types';
import { getWorkflowBlock, isTopLevelBlock } from '../workflow';

// What a phase may need to know about the run, beyond the desc itself. A phase
// reads only the parts it cares about (most ignore it).
export interface WorkflowMutatorContext {
  // The project the document lives in.
  projectName: string;
}

// Base for one phase of normalizing a workflow document after an in-memory
// change. A subclass implements mutate() for its phase; the base offers block
// read/write helpers so a phase works at the desc level while the optional
// { workflow } nesting is handled in one place.
export abstract class WorkflowMutatorBase {
  // A label for the phase (useful for logging / ordering).
  abstract readonly name: string;

  // Runs this phase, returning the corrected desc (or the same one unchanged).
  abstract mutate(desc: WorkflowDesc, context: WorkflowMutatorContext): WorkflowDesc;

  // The workflow block within a desc (unwraps the optional { workflow } level).
  protected getBlock(desc: WorkflowDesc): WorkflowBlock | undefined {
    return getWorkflowBlock(desc.workflow);
  }

  // Writes a block back into a desc, preserving the nesting the desc used.
  protected writeBlock(desc: WorkflowDesc, block: WorkflowBlock): WorkflowDesc {
    if (isTopLevelBlock(desc.workflow)) {
      return { ...desc, workflow: block };
    }
    return { ...desc, workflow: { ...desc.workflow, workflow: block } };
  }
}
