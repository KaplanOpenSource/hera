import { DocumentObj } from '../objects/DocumentObj';
import { ProjectObj } from '../objects/ProjectObj';
import { VALUE_GROUP_UNDEFINED } from '../utils/splitTree';
import { CENTRAL_REPO_FOLDER_ID, idFromDocId, idFromRepoId, isSplitId, toolkitNameFromSplitId } from './idDocId';
import { isWorkflowDoc } from './workflow';

export enum TabKind {
  Notebook = 'notebook',
  Document = 'document',
  Agent = 'agent',
  Workflow = 'workflow',
  ProjectConfig = 'projectConfig',
  Repository = 'repository',
  CentralRepository = 'centralRepository',
  Toolkit = 'toolkit',
}

const isAgentDoc = (doc: DocumentObj) => {
  const resource = doc.data.resource;
  return typeof resource === 'object' && resource.effects !== undefined;
};

// Classifies a document straight from its DocumentObj (the tree already holds one),
// without needing an item id / project lookup like classifyTab does.
export const classifyDocument = (doc: DocumentObj): TabKind => {
  if (doc.isConfig) return TabKind.ProjectConfig;
  if (doc.isNotebook) return TabKind.Notebook;
  if (isWorkflowDoc(doc.data)) return TabKind.Workflow;
  if (isAgentDoc(doc)) return TabKind.Agent;
  return TabKind.Document;
};

export const classifyTab = (showItemId: string, project: ProjectObj): TabKind | undefined => {
  if (showItemId === CENTRAL_REPO_FOLDER_ID) return TabKind.CentralRepository;

  if (isSplitId(showItemId)) return TabKind.Toolkit;

  if (idFromRepoId(showItemId)) return TabKind.Repository;

  const docid = idFromDocId(showItemId);
  if (docid) {
    const doc = project.allDocuments.find(d => d.docid === docid);
    if (!doc) return undefined;
    if (doc.isConfig) return TabKind.ProjectConfig;
    if (doc.isNotebook) return TabKind.Notebook;
    if (isWorkflowDoc(doc.data)) return TabKind.Workflow;
    if (isAgentDoc(doc)) return TabKind.Agent;
    return TabKind.Document;
  }

  return TabKind.ProjectConfig;
};

export const tabKindClassName = (showItemId: string, project: ProjectObj): string | undefined => {
  const kind = classifyTab(showItemId, project);
  return kind ? `tab-kind-${kind}` : undefined;
};

// Returns a human-readable tab title for a given item ID.
export const detailsTabName = (showItemId: string, project: ProjectObj): string => {
  if (showItemId === CENTRAL_REPO_FOLDER_ID) return 'Repositories';
  if (isSplitId(showItemId)) {
    const toolkitName = toolkitNameFromSplitId(showItemId, project.documents);
    if (toolkitName) return toolkitName === VALUE_GROUP_UNDEFINED ? 'No toolkit' : toolkitName;
  }
  const docid = idFromDocId(showItemId);
  if (docid) {
    const doc = project.allDocuments.find(d => d.docid === docid);
    if (doc) return doc.isConfig ? project.name + ' config' : doc.name;
  }
  const repoPath = idFromRepoId(showItemId);
  if (repoPath) return repoPath;
  return project.name + ' config';
};
