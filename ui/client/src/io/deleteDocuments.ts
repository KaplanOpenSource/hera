import { ProjectEntire } from "@shared/types";
import { useProjectStore } from "../stores/useProjectStore";
import { postJson } from "./entrypoint";

// Deletes the given documents. The answer carries the project that is left, so
// nothing is reloaded after it. Returns false on failure.
export const deleteDocuments = async (docOids: string[]) => {
  const data = await postJson('/document/delete', { docIds: docOids },
    { label: `delete ${docOids.length} documents` });
  if (!data) {
    return false;
  }
  useProjectStore.getState().setCurrentProject(data.project as ProjectEntire);
  return true;
};
