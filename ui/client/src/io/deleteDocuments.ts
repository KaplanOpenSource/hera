import { useProjectStore } from "../stores/useProjectStore";
import { documentApi } from "../api/api";

// Deletes the given documents. The answer carries the project that is left, so
// nothing is reloaded after it. Returns false on failure.
export const deleteDocuments = async (docOids: string[]) => {
  const data = await documentApi.remove(docOids);
  if (!data) {
    return false;
  }
  useProjectStore.getState().setCurrentProject(data.project);
  return true;
};
