import { DocumentObj } from '../objects/DocumentObj';
import { ProjectObj } from '../objects/ProjectObj';
import { CENTRAL_REPO_FOLDER_ID, idFromDocId, idFromRepoId } from './idDocId';

export enum TabKind {
  Notebook = 'notebook',
  Document = 'document',
  Agent = 'agent',
  ProjectConfig = 'projectConfig',
  Repository = 'repository',
  CentralRepository = 'centralRepository',
}

const isAgentDoc = (doc: DocumentObj) => {
  const resource = doc.data.resource;
  return typeof resource === 'object' && resource.effects !== undefined;
};

export const classifyTab = (showItemId: string, project: ProjectObj): TabKind | undefined => {
  if (showItemId === CENTRAL_REPO_FOLDER_ID) return TabKind.CentralRepository;

  if (idFromRepoId(showItemId)) return TabKind.Repository;

  const docid = idFromDocId(showItemId);
  if (docid) {
    const doc = project.allDocuments.find(d => d.docid === docid);
    if (!doc) return undefined;
    if (doc.isConfig) return TabKind.ProjectConfig;
    if (doc.isNotebook) return TabKind.Notebook;
    if (isAgentDoc(doc)) return TabKind.Agent;
    return TabKind.Document;
  }

  return TabKind.ProjectConfig;
};
