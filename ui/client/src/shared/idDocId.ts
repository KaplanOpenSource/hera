import { getValueAtPath } from '../utils/compareJsons';
import { VALUE_GROUP_UNDEFINED } from '../utils/splitTree';

const ID_PREFIX_DOC = 'document_'
const ID_PREFIX_REPO = '__repo_*_'
const ID_PREFIX_SPLIT = 'split_'
const ID_PREFIX_TOOLKIT = 'split_/toolkit='
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

export const idFromToolkitSplitId = (id: string) => {
  return id && id.startsWith(ID_PREFIX_TOOLKIT) ? id.replace(ID_PREFIX_TOOLKIT, "") : undefined;
};

export const isSplitId = (id: string) => {
  return id?.startsWith(ID_PREFIX_SPLIT);
};

// Resolves the toolkit name for any split tree node ID.
// For toolkit splits (split_/toolkit=X), returns X directly.
// For child splits (split_/type=Y), finds a matching document and returns its toolkit.
export const toolkitNameFromSplitId = (
  splitId: string,
  documents: { extDesc: Record<string, any>, toolkit: string | undefined }[],
): string | undefined => {
  const toolkitName = idFromToolkitSplitId(splitId);
  if (toolkitName) return toolkitName;

  const match = splitId.match(/^split_(.+)=(.+)$/);
  if (!match) return undefined;
  const [, path, value] = match;
  const doc = documents.find(d =>
    String(getValueAtPath(d.extDesc, path)) === value
  );
  return doc?.toolkit ?? VALUE_GROUP_UNDEFINED;
};
