const ID_PREFIX_DOC = 'document_'
const ID_PREFIX_REPO = '__repo_*_'
export const ID_NOTEBOOKS_GROUP = 'notebooks-group'
const ID_PREFIX_NOTEBOOK = 'notebook_'
export const TEMP_REPO_NAME = '*Temp Repository*';

export const idDocId = (oid: string) => {
  return ID_PREFIX_DOC + oid;
}

export const idFromDocId = (docid: string) => {
  return docid && docid.startsWith(ID_PREFIX_DOC) ? docid.replace(ID_PREFIX_DOC, "") : undefined;
}

export const idRepoId = (repoName: string) => {
  return ID_PREFIX_REPO + repoName;
}

export const idFromRepoId = (repoid: string) => {
  return repoid && repoid.startsWith(ID_PREFIX_REPO) ? repoid.replace(ID_PREFIX_REPO, "") : undefined;
};

export const idNotebookId = (name: string) => {
  return ID_PREFIX_NOTEBOOK + name;
}

export const idFromNotebookId = (itemId: string) => {
  return itemId && itemId.startsWith(ID_PREFIX_NOTEBOOK) ? itemId.replace(ID_PREFIX_NOTEBOOK, "") : undefined;
};
