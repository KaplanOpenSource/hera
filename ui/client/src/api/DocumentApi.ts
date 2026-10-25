import { ApiCall } from './ApiCall';
import { Collection, DocumentDesc, ProjectDocument, ProjectEntire } from '../shared/types';
import { FORBIDDEN_FIELDS } from '../shared/constants';

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

export type UpdatedDocument = {
  document: ProjectDocument,
  project: ProjectEntire,
};

// Every call to /document/*. Each method answers data, or null when the call
// failed.
export class DocumentApi extends ApiCall {
  // Adds one document and answers it with the project it joined, or null when
  // the call failed.
  async create(newDocument: NewDocument, label: string): Promise<CreatedDocument | null> {
    const data = await this.postJson('/document/create', newDocument, { label });
    if (!data) {
      return null;
    }
    return { project: data.project as ProjectEntire, addedDocId: data.addedDocId as string };
  }

  // Writes the whole document, minus the fields that say which document it is.
  // The answer carries the project too, so saving is one call. Returns null
  // when the call failed.
  async update(newDoc: ProjectDocument): Promise<UpdatedDocument | null> {
    const docId = newDoc._id.$oid;
    const fields: { [field: string]: any } = {};
    for (const [field, value] of Object.entries(newDoc)) {
      if (!FORBIDDEN_FIELDS.includes(field)) {
        fields[field] = value;
      }
    }
    const data = await this.postJson('/document/update', { docId, fields },
      { label: 'update document' });
    if (!data) {
      return null;
    }
    return { document: data.document as ProjectDocument, project: data.project as ProjectEntire };
  }

  // Deletes the given documents. The answer carries the project that is left.
  // Returns null when the call failed.
  async remove(docOids: string[]): Promise<{ project: ProjectEntire } | null> {
    const data = await this.postJson('/document/delete', { docIds: docOids },
      { label: `delete ${docOids.length} documents` });
    if (!data) {
      return null;
    }
    return { project: data.project as ProjectEntire };
  }
}
