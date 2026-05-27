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
  if (!id?.startsWith(ID_PREFIX_TOOLKIT)) return undefined;
  const rest = id.slice(ID_PREFIX_TOOLKIT.length);
  const slashIdx = rest.indexOf('/');
  return slashIdx === -1 ? rest : rest.slice(0, slashIdx);
};

export const isSplitId = (id: string) => {
  return id?.startsWith(ID_PREFIX_SPLIT);
};

// Resolves the toolkit name for any split tree node ID.
// For toolkit splits (split_/toolkit=X), returns X directly.
// For child splits (split_/toolkit=X/...), extracts X from the key prefix.
// For splits without a toolkit ancestor, finds the toolkit from matching documents.
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

// Normalizes any split ID to its parent toolkit split ID.
// Child splits reuse the parent toolkit's details tab instead of opening a new one.
export const normalizeSplitId = (
  splitId: string,
  documents: { extDesc: Record<string, any>, toolkit: string | undefined }[],
): string => {
  const toolkitName = toolkitNameFromSplitId(splitId, documents);
  return ID_PREFIX_TOOLKIT + (toolkitName ?? VALUE_GROUP_UNDEFINED);
};
