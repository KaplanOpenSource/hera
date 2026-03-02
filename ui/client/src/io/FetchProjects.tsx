import { ProjectEntire, ProjectName, Toolkit } from "@shared/types";
import { useEffect } from "react";
import { useNavigate } from "react-router-dom";
import { DEFAULT_PROJECT, NO_PROJECT, useProjectStore } from "../stores/useProjectStore";
import { execPython } from "./execPython";
import { fetchPython } from "./fetchPython";
import { ProjectCommands } from "./ProjectCommands";

export const fetchProjectsNames = async () => {
  const { data, problem } = await fetchPython(ProjectCommands.projectNames());

  if (!problem) {
    const projects = (data.projects || []) as ProjectName[];
    const first = projects.filter(({ name }) => name === DEFAULT_PROJECT);
    const rest = projects.filter(({ name }) => name !== DEFAULT_PROJECT);
    rest.sort((a, b) => a.name.localeCompare(b.name, undefined, { sensitivity: 'base' }));
    useProjectStore.getState().setProjectNames([...first, ...rest])
  }
}

export const fetchProjectDetails = async (projectName: string) => {
  const { setCurrentProject } = useProjectStore.getState();
  const { data, problem } = await execPython(`
import json
from hera.datalayer import All
docs = All.getDocumentsAsDict('${projectName}', with_id=True)
project = {"name": '${projectName}', "documents": docs['documents']}
result = json.dumps(project,indent=4)
  `);
  if (!problem) {
    const project = JSON.parse(data) as ProjectEntire;
    setCurrentProject(project);
  }
}

export const fetchToolkits = async (projectName: string) => {
  const { setToolkits: setToolKits } = useProjectStore.getState();
  const { data, problem } = await execPython(`
from hera import toolkitHome
import json
result = toolkitHome.getToolkitDocuments()
`);
  // table = toolkitHome.getToolkitTable('${name}')
  // result = table.to_json(orient='records')
  if (!problem) {
    const newToolkits = data.map((d: any) => {
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
    // setToolKits(JSON.parse(data) as Toolkit[])
    setToolKits(newToolkits)
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
  }, [])

  // Sync project selection with URL and fetch project data
  useEffect(() => {
    if (projectNames.length === 0) return;

    const urlProject = resolveProjectFromUrl(urlProjectName, projectNames);

    // URL changed to a different valid project → switch to it
    if (urlProject && urlProject !== currProjectName) {
      selectProject(urlProject);
      return;
    }

    if (currProjectName === NO_PROJECT) {
      selectProject(urlProject ?? projectNames[0].name);
    } else {
      fetchToolkits(currProjectName);
      fetchProjectDetails(currProjectName);
      const expectedPath = '/' + encodeURIComponent(currProjectName);
      if (location.pathname !== expectedPath) {
        navigate(expectedPath, { replace: true });
      }
    }
  }, [currProjectName, projectNames, urlProjectName, navigate]);

  return null;
}

