import { useProjectStore } from "../stores/useProjectStore";
import { fetchProjectDetails } from "./FetchProjects";
import { fetchPython } from "./fetchPython";

// Deletes the given documents, then reloads the project. Returns false on failure.
export const deleteDocuments = async (docOids: string[]) => {
  const projectName = useProjectStore.getState().currProjectName;
  const { data } = await fetchPython({
    results: [],
    label: `delete ${docOids.length} documents`,
    code: [
      'from hera.datalayer import All',
      ...docOids.map(oid => `All.deleteDocumentByID('${oid}')`),
    ].join('\n'),
  });
  if (!data) {
    return false;
  }
  await fetchProjectDetails(projectName);
  return true;
};
