import { ProjectEntire } from '@shared/types';
import { projectApi } from '../api/api';
import { ProjectObj } from '../objects/ProjectObj';
import { useProjectStore } from '../stores/useProjectStore';
import { pushInfo } from './snackbar';

// Registers the untracked notebooks in the project's files directory; runs after a project loads.
export const detectNotebooks = async (project: ProjectEntire) => {
  const filesDir = new ProjectObj(project).configDocument?.data.desc.filesDirectory;
  if (!filesDir) {
    return;
  }
  const data = await projectApi.detectNotebooks(project.name, filesDir);
  if (!data) {
    return;
  }
  const added = data.addedNotebooks as number;
  if (added === 0) {
    return;
  }
  const { currProjectName, setCurrentProject } = useProjectStore.getState();
  if (currProjectName !== project.name) {
    return;
  }
  setCurrentProject(data.project as ProjectEntire);
  let plural = 's';
  if (added === 1) {
    plural = '';
  }
  pushInfo(`Detected ${added} new notebook${plural}`);
};
