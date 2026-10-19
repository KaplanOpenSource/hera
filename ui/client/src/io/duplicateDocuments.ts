import { Collection, ProjectDocument } from "../shared/types";
import { useProjectStore } from "../stores/useProjectStore";
import { documentApi } from "../api/api";
import { NewDocument } from "../api/DocumentApi";

// The document's _cls says which collection it lives in, e.g. Metadata.Cache.
const COLLECTION_OF_CLS: { [cls: string]: Collection } = {
  Measurements: Collection.Measurements,
  Simulations: Collection.Simulations,
  Cache: Collection.Cache,
};

const collectionOfCls = (cls: string) => {
  return COLLECTION_OF_CLS[cls?.split('.').pop() ?? ''];
};

// The first "<name> copy", "<name> copy 2", ... that no document in the project uses.
export const freeCopyName = (name: string, taken: string[]) => {
  let copyName = `${name} copy`;
  let num = 2;
  while (taken.includes(copyName)) {
    copyName = `${name} copy ${num}`;
    num += 1;
  }
  return copyName;
};

const copyOf = (doc: ProjectDocument, taken: string[]): NewDocument | null => {
  const collection = collectionOfCls(doc._cls);
  if (!collection) {
    console.error(`cannot duplicate a ${doc._cls} document`);
    return null;
  }
  const desc = { ...doc.desc };
  if (desc.datasourceName) {
    desc.datasourceName = freeCopyName(desc.datasourceName, taken);
  }
  return {
    projectName: doc.projectName,
    collection,
    type: doc.type,
    resource: doc.resource,
    dataFormat: doc.dataFormat,
    desc,
  };
};

// Copies each document into its own collection, with a free "copy" name. The
// copies come from the store, so nothing is read back first. Returns false on
// failure.
export const duplicateDocuments = async (docOids: string[]) => {
  const documents = useProjectStore.getState().currProject?.documents ?? [];
  const taken = documents.map(d => d.desc?.datasourceName).filter(Boolean) as string[];
  for (const oid of docOids) {
    const doc = documents.find(d => d._id.$oid === oid);
    if (!doc) {
      console.error(`no document ${oid} in the project`);
      return false;
    }
    const copy = copyOf(doc, taken);
    if (!copy) {
      return false;
    }
    // Each copy's name is taken by the next one, so two copies never collide.
    const copyName = copy.desc?.datasourceName;
    if (copyName) {
      taken.push(copyName);
    }
    const created = await documentApi.create(copy, `duplicate document ${doc.desc?.datasourceName ?? oid}`);
    if (!created) {
      return false;
    }
    // The last answer is the freshest project, so storing each is harmless.
    useProjectStore.getState().setCurrentProject(created.project);
  }
  return true;
};
