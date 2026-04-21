import { ProjectEntire, ProjectName, Toolkit } from "@shared/types";
import { useEffect } from "react";
import { useNavigate } from "react-router-dom";
import { NO_PROJECT, useProjectStore } from "../stores/useProjectStore";
import { fetchPython } from "./fetchPython";
import { ProjectCommands } from "./ProjectCommands";

export const fetchProjectsNames = async () => {
  const { data } = await fetchPython(ProjectCommands.projectNames());
  if (data) {
    useProjectStore.getState().setProjectNames(data.projects as ProjectName[])
  }
}

export const fetchProjectDetails = async (projectName: string) => {
  const { data } = await fetchPython({
    results: ['project'],
    label: `project details ${projectName}`,
    code: `
from hera.datalayer import All
docs = All.getDocumentsAsDict('${projectName}', with_id=True)
project = {"name": '${projectName}', "documents": docs['documents']}
`,
  });
  if (data) {
    // Skip if the user switched projects while this fetch was in flight
    const { currProjectName, setCurrentProject } = useProjectStore.getState();
    if (currProjectName === projectName) {
      setCurrentProject(data.project as ProjectEntire);
    }
  }
}

const parseToolkits = (toolkitDocs: any[]): Toolkit[] =>
  toolkitDocs.map((d: any) => {
    const desc = d.desc;
    const t = {
      toolkit: d.toolkit,
      cls: desc.classpath,
    } as Toolkit;
    for (const field of [
      'description',
      'source',
      'type',
      'repositoryName',
      'version',
    ] as const) {
      if (desc[field]) t[field] = desc[field];
    }
    return t;
  });

export const fetchProjectData = async (projectName: string) => {
  const { data } = await fetchPython(
    {
      results: ['toolkitDocs'],
      label: 'toolkits',
      code: `
from hera import toolkitHome
toolkitDocs = toolkitHome.getToolkitDocuments()
`,
    },
    {
      results: ['project'],
      label: `project ${projectName}`,
      code: `
from hera.datalayer import All
docs = All.getDocumentsAsDict('${projectName}', with_id=True)
project = {"name": '${projectName}', "documents": docs['documents']}
`,
    },
  );
  if (data) {
    const { currProjectName, setCurrentProject, setToolkits } = useProjectStore.getState();
    setToolkits(parseToolkits(data.toolkitDocs));
    if (currProjectName === projectName) {
      setCurrentProject(data.project as ProjectEntire);
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
  const {
    projectNames,
    currProjectName,
    selectProject,
  } = useProjectStore();
  const navigate = useNavigate();

  useEffect(() => {
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
      fetchProjectData(currProjectName);
      const expectedPrefix = '/' + encodeURIComponent(currProjectName);
      if (!location.pathname.startsWith(expectedPrefix)) {
        navigate(expectedPrefix, { replace: true });
      }
    }
  }, [currProjectName, projectNames, urlProjectName, navigate]);

  return null;
}

