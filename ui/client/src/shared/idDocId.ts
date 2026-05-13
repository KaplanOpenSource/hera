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
