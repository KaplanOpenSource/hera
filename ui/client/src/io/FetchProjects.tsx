import { ProjectEntire, ProjectName } from "@shared/types";
import { useEffect, useRef } from "react";
import { useNavigate } from "react-router-dom";
import { NO_PROJECT, useProjectStore } from "../stores/useProjectStore";
import { useProjectListStore } from "../stores/useProjectListStore";
import { useToolkitStore } from "../stores/useToolkitStore";
import { getJson } from "./entrypoint";
import { parseToolkits } from "./parseToolkits";

export const fetchProjectsNames = async () => {
  const data = await getJson('/project/names', {}, { label: 'project names' });
  if (data) {
    useProjectListStore.getState().setProjectNames(data.projects as ProjectName[])
  }
}

// silent: no toasts (e.g. for the periodic auto-reload).
export const fetchProjectDetails = async (projectName: string, silent = false) => {
  const data = await getJson('/project/documents', { projectName },
    { label: `project details ${projectName}`, silent });
  if (data) {
    // Skip if the user switched projects while this fetch was in flight
    const { currProjectName, setCurrentProject } = useProjectStore.getState();
    if (currProjectName === projectName) {
      setCurrentProject(data.project as ProjectEntire);
    }
  }
}

export const fetchProjectData = async (projectName: string) => {
  const [toolkits, details] = await Promise.all([
    getJson('/project/toolkits', {}, { label: 'toolkits' }),
    getJson('/project/documents', { projectName }, { label: `project ${projectName}` }),
  ]);
  if (toolkits) {
    useToolkitStore.getState().setToolkits(parseToolkits(toolkits.toolkitDocs));
  }
  if (details) {
    const { currProjectName, setCurrentProject } = useProjectStore.getState();
    if (currProjectName === projectName) {
      setCurrentProject(details.project as ProjectEntire);
    }
  }
}

export const resolveProjectFromUrl = (urlProjectName: string | undefined, projectNames: ProjectName[]): string | undefined => {
  if (!urlProjectName) return undefined;
  const decodedName = decodeURIComponent(urlProjectName);
  return projectNames.find(p => p.name === decodedName)?.name;
};

export const FetchProjects = ({
  urlProjectName,
}: {
  urlProjectName?: string,
}) => {
  const { projectNames } = useProjectListStore();
  const { currProjectName, selectProject } = useProjectStore();
  const navigate = useNavigate();
  const firstRun = useRef(true);

  // Startup brought the names; refresh them when the URL project changes, so a
  // project just added is in the list.
  useEffect(() => {
    if (firstRun.current) {
      firstRun.current = false;
      return;
    }
    fetchProjectsNames();
  }, [urlProjectName])

  // Sync project selection with URL and fetch project data
  useEffect(() => {
    if (projectNames.length === 0) return;

    const urlProject = resolveProjectFromUrl(urlProjectName, projectNames);

    // URL changed to a different valid project → switch to it
    if (urlProject && urlProject !== currProjectName) {
      selectProject(urlProject);
      return;
    }

    // URL points to a project not yet in the list — wait for list to update
    if (urlProjectName && !urlProject) {
      return;
    }

    if (currProjectName === NO_PROJECT) {
      selectProject(urlProject ?? projectNames[0].name);
    } else {
      // Skip the fetch when the store already holds this project, e.g. startup
      // loaded it (see fetchInit).
      if (useProjectStore.getState().currProject?.name !== currProjectName) {
        fetchProjectData(currProjectName);
      }
      const expectedPrefix = '/' + encodeURIComponent(currProjectName);
      if (!location.pathname.startsWith(expectedPrefix)) {
        navigate(expectedPrefix, { replace: true });
      }
    }
  }, [currProjectName, projectNames, urlProjectName, navigate]);

  return null;
}
