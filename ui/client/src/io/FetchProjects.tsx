import { ProjectEntire, ProjectName, Toolkit } from "@shared/types";
import { DEFAULT_PROJECT, NO_PROJECT, useProjectStore } from "../stores/useProjectStore";
import { execPython } from "./execPython";
import { useEffect } from "react";

export const fetchProjectsNames = async () => {
  const { data, problem } = await execPython(`
from hera.datalayer.project import getProjectList
result = [{"name": proj} for proj in getProjectList()]
  `);

  if (!problem) {
    const projects = (data || []) as ProjectName[];
    const first = projects.filter(({ name }) => name === DEFAULT_PROJECT);
    const rest = projects.filter(({ name }) => name !== DEFAULT_PROJECT);
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
table = toolkitHome.getToolkitTable('${projectName}')
result = table.to_json(orient='records', indent=2)
  `);
  if (!problem) {
    setToolKits(JSON.parse(data) as Toolkit[])
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

