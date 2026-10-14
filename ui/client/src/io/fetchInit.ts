import { ProjectEntire, ProjectName } from '@shared/types';
import { BASEURL } from '../shared/baseurl';
import { useProjectListStore } from '../stores/useProjectListStore';
import { useProjectStore } from '../stores/useProjectStore';
import { useServerConstants } from '../stores/useServerConstants';
import { useServerUserStore } from '../stores/useServerUserStore';
import { useToolkitStore } from '../stores/useToolkitStore';
import { parseToolkits } from './parseToolkits';

// What GET /init answers. Any piece may be missing: the server collects each
// one on its own and names the failures in `problems` (see init_data.py).
export type InitData = {
  ready: boolean,
  projects?: ProjectName[],
  datatypes?: { [key: string]: string },
  toolkitDocs?: any[],
  project?: ProjectEntire,
  server?: { heraVersion?: string, username?: string, inDocker?: boolean },
  problems?: { [piece: string]: string },
};

// Asks the server for everything the UI needs to start. Returns the answer, or
// null when the server did not answer at all. Not ready yet carries no data.
export const fetchInit = async (wantedProjectName?: string): Promise<InitData | null> => {
  const query = wantedProjectName ? `?project=${encodeURIComponent(wantedProjectName)}` : '';
  const response = await fetch(`${BASEURL}/init${query}`);
  const contentType = response.headers.get('content-type');
  // Old server (no /init route) answers with the SPA's HTML — treat as no answer.
  if (!response.ok || !contentType?.includes('application/json')) {
    return null;
  }
  const data = await response.json() as InitData;
  if (!data.ready) {
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
  }
  return data;
};
