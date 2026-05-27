const ID_PREFIX_DOC = 'document_'
const ID_PREFIX_REPO = '__repo_*_'
export const TEMP_REPO_NAME = '*Temp Repository*';
export const CENTRAL_REPO_FOLDER_ID = 'central-repo-folder';

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

const ID_PREFIX_TOOLKIT = 'split_/toolkit='

export const idFromToolkitSplitId = (id: string) => {
  return id && id.startsWith(ID_PREFIX_TOOLKIT) ? id.replace(ID_PREFIX_TOOLKIT, "") : undefined;
};
