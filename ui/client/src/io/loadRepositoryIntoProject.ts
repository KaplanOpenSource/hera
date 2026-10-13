import { useProjectStore } from "../stores/useProjectStore";
import { fetchProjectDetails } from "./FetchProjects";
import { fetchPython } from "./fetchPython";

// Loads all the datasources of a repository into the current project.
export const loadRepositoryIntoProject = async (repositoryName: string) => {
  const projectName = useProjectStore.getState().currProjectName;
  await fetchPython({
    results: [],
    label: `load repository ${repositoryName}`,
    code: `
from hera.utils.data.toolkit import dataToolkit
dataToolkit().loadAllDatasourcesInRepositoryToProject(projectName='${projectName}', repositoryName='${repositoryName}', overwrite=False)
`,
  });
  await fetchProjectDetails(projectName);
};
