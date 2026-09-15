import { NO_PROJECT, useProjectStore } from '../../stores/useProjectStore';
import { WorkflowDesc } from '../types';
import { fillProjectName } from './mutators/FillProjectNameMutator';
import { syncParameters } from './mutators/SyncParametersMutator';

export class MutatorsListHandler {
  // Corrects a workflow desc after an in-memory change (edit and create, not
  // load): seed every empty project-name field with the current project, then
  // rebuild the parameters index from the result. A non-workflow desc comes
  // back with no index, since it has no workflow block.
  static normalize(desc: WorkflowDesc): WorkflowDesc {
    let result = desc;
    const projectName = useProjectStore.getState().currProjectName;
    if (projectName && projectName !== NO_PROJECT) {
      result = fillProjectName(result, projectName);
    }
    return syncParameters(result);
  }
}
