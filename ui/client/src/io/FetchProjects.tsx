import { ProjectEntire, ProjectName, Toolkit } from "@shared/types";
import { useEffect } from "react";
import { DEFAULT_PROJECT, NO_PROJECT, useProjectStore } from "../stores/useProjectStore";
import { execPython } from "./execPython";

export const fetchProjectsNames = async () => {
  const { data, problem } = await execPython(`
from hera.datalayer.project import getProjectList
result = [{"name": proj} for proj in getProjectList()]
  `);

  if (!problem) {
    const projects = (data || []) as ProjectName[];
    const first = projects.filter(({ name }) => name === DEFAULT_PROJECT);
    const rest = projects.filter(({ name }) => name !== DEFAULT_PROJECT);
    rest.sort((a, b) => a.name.localeCompare(b.name, undefined, { sensitivity: 'base' }));
    useProjectStore.getState().setProjectNames([...first, ...rest])
    // console.log('inside:', useProjectStore.getState().projectNames, data)
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

export const FetchProjects = ({ }) => {
  const {
    projectNames,
    currProjectName,
    selectProject,
  } = useProjectStore();

  useEffect(() => {
    fetchProjectsNames();
  }, [])

  // console.log('outside:', useProjectStore.getState().projectNames)

  useEffect(() => {
    if (currProjectName === NO_PROJECT && projectNames.length > 0) {
      selectProject(projectNames[0].name);
    } else if (currProjectName !== NO_PROJECT) {
      fetchToolkits(currProjectName);
      fetchProjectDetails(currProjectName);
    }
  }, [currProjectName, projectNames]);

  return null;
}

