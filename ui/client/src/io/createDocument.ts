import { Collection, DocumentDesc, ProjectEntire } from "../shared/types";
import { postJson } from "./entrypoint";

// A notebook document owns an .ipynb file; the server writes it if it is
// missing, with one toolkit line per name.
export type NewNotebook = {
  toolkitNames: string[],
};

export type NewDocument = {
  projectName: string,
  collection: Collection,
  type?: string,
  // A path, or the data itself (an agent's resource is an object).
  resource?: unknown,
  dataFormat?: string,
  desc?: DocumentDesc,
  notebook?: NewNotebook,
};

export type CreatedDocument = {
  project: ProjectEntire,
  addedDocId: string,
};

// Adds one document and answers it with the project it joined, or null when the
// call failed.
export const createDocument = async (
  newDocument: NewDocument,
  label: string,
): Promise<CreatedDocument | null> => {
  const data = await postJson('/document/create', newDocument, { label });
  if (!data) {
    return null;
  }
  return { project: data.project as ProjectEntire, addedDocId: data.addedDocId as string };
};
