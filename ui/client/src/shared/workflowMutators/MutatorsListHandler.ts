import { WorkflowDesc } from '../types';
import { FillProjectNameMutator } from './mutators/FillProjectNameMutator';
import { SyncParametersMutator } from './mutators/SyncParametersMutator';
import { WorkflowMutatorBase } from './WorkflowMutatorBase';

// Holds the workflow mutators and runs them as one pipeline over a workflow
// desc, on every in-memory change (edit and create, not load).
export class MutatorsListHandler {
  // The phases, in run order. Add a mutator here to include it in the pipeline.
  static readonly mutators: WorkflowMutatorBase[] = [
    new FillProjectNameMutator(),
    new SyncParametersMutator(),
  ];

  // Runs every mutator in turn over the desc. A non-workflow desc is returned
  // unchanged (each phase no-ops without a workflow block).
  static normalize(desc: WorkflowDesc, projectName: string): WorkflowDesc {
    let result = desc;
    for (const mutator of MutatorsListHandler.mutators) {
      result = mutator.mutate(result, { projectName });
    }
    return result;
  }
}
