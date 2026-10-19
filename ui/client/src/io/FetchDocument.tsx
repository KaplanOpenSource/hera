import { ProjectDocument, ProjectEntire } from "@shared/types";
import { postJson } from "./entrypoint";
import { FORBIDDEN_FIELDS } from "../shared/constants";

export type UpdatedDocument = {
  document: ProjectDocument,
  project: ProjectEntire,
};

// Writes the whole document, minus the fields that say which document it is.
// The answer carries the project too, so saving is one call. Returns null when
// the call failed.
export const updateDocument = async (newDoc: ProjectDocument): Promise<UpdatedDocument | null> => {
  const docId = newDoc._id.$oid;
  const fields: { [field: string]: any } = {};
  for (const [field, value] of Object.entries(newDoc)) {
    if (!FORBIDDEN_FIELDS.includes(field)) {
      fields[field] = value;
    }
  }
  const data = await postJson('/document/update', { docId, fields },
    { label: 'update document' });
  if (!data) {
    return null;
  }
  return { document: data.document as ProjectDocument, project: data.project as ProjectEntire };
}
