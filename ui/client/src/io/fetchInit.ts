import { serverApi } from '../api/api';
import { InitData } from '../api/ServerApi';
import { useProjectListStore } from '../stores/useProjectListStore';
import { useProjectStore } from '../stores/useProjectStore';
import { useServerConstants } from '../stores/useServerConstants';
import { useServerUserStore } from '../stores/useServerUserStore';
import { useToolkitStore } from '../stores/useToolkitStore';
import { parseToolkits } from '../utils/parseToolkits';
import { detectNotebooks } from './detectNotebooks';

// Asks the server for everything the UI needs to start. Returns the answer, or
// null when the server did not answer at all. Not ready yet carries no data.
export const fetchInit = async (wantedProjectName?: string): Promise<InitData | null> => {
  const data = await serverApi.init(wantedProjectName);
  if (!data || !data.ready) {
    return data;
  }
  if (data.projects) {
    useProjectListStore.getState().setProjectNames(data.projects);
  }
  if (data.datatypes) {
    useServerConstants.getState().setDataTypes(data.datatypes);
  }
  if (data.toolkitDocs) {
    useToolkitStore.getState().setToolkits(parseToolkits(data.toolkitDocs));
  }
  if (data.server) {
    useServerUserStore.getState().setServerInfo(data.server);
  }
  if (data.project) {
    const { selectProject, setCurrentProject } = useProjectStore.getState();
    selectProject(data.project.name);
    setCurrentProject(data.project);
    // Not awaited, so the project shows before the scan ends.
    detectNotebooks(data.project);
  }
  return data;
};
