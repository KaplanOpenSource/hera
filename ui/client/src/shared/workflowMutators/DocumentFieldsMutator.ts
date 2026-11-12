import { NO_PROJECT, useProjectStore } from '../../stores/useProjectStore';
import { WorkflowDesc } from '../types';
import { fillProjectName } from './mutators/FillProjectNameMutator';
import { syncParameters } from './mutators/SyncParametersMutator';

export class DocumentFieldsMutator {
  // Corrects a workflow desc after an in-memory change: fill every empty
  // project-name field, then rebuild the parameters index from the result.
  static mutate(desc: WorkflowDesc): WorkflowDesc {
    let result = desc;
    const projectName = useProjectStore.getState().currProjectName;
    if (projectName && projectName !== NO_PROJECT) {
      result = fillProjectName(result, projectName);
    }
    return syncParameters(result);
  }
}
