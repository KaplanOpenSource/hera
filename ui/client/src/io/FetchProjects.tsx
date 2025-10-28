import { ProjectEntire, ProjectName } from "@shared/types";
import { DEFAULT_PROJECT, NO_PROJECT, useProjectStore } from "../stores/useProjectStore";
import { execPython } from "./execPython";
import { useEffect } from "react";

export const fetchProjectsNames = async () => {
  const { setProjectNames } = useProjectStore.getState();
  const { data, problem } = await execPython(`
from hera.datalayer.project import getProjectList
result = [{"name": proj} for proj in getProjectList()]
  `);

  if (!problem) {
    const projects = (data || []) as ProjectName[];
    const first = projects.filter(({ name }) => name === DEFAULT_PROJECT);
    const rest = projects.filter(({ name }) => name !== DEFAULT_PROJECT);
    setProjectNames([...first, ...rest])
  }
}

export const fetchProjectDetails = async (projectName: string) => {
  const { setCurrentProject } = useProjectStore.getState();
  const { data, problem } = await execPython(`
import json
from hera.datalayer import All
docList = []
for doc in All.getDocuments('${projectName}'):
    docDict = doc.asDict()
    if ('docid' not in docDict['desc']):
        docDict['desc']['docid'] = str(doc.id)
    docList.append(docDict)
project = {"name": '${projectName}', "documents": docList}
result = json.dumps(project,indent=4)
  `);
  if (!problem) {
    const project = JSON.parse(data) as ProjectEntire;
    setCurrentProject(project);
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

  useEffect(() => {
    if (currProjectName === NO_PROJECT && projectNames.length > 0) {
      selectProject(projectNames[0].name);
    } else if (currProjectName !== NO_PROJECT) {
      fetchProjectDetails(currProjectName);
    }
  }, [currProjectName, projectNames]);

  return null;
}

