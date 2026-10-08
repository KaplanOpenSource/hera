import { useProjectStore } from "../stores/useProjectStore";
import { fetchProjectDetails } from "./FetchProjects";
import { fetchPython } from "./fetchPython";

// Copies a document into the same collection, with a free "copy" name.
const DUPLICATE_HELPER = `
from hera.datalayer import All

def duplicateDocumentByID(oid):
    doc = All.getDocumentByID(oid)
    desc = dict(doc.desc)
    name = desc.get('datasourceName')
    if name:
        taken = [d.desc.get('datasourceName') for d in All.getDocuments(projectName=doc.projectName)]
        newName = name + ' copy'
        num = 2
        while newName in taken:
            newName = name + ' copy ' + str(num)
            num += 1
        desc['datasourceName'] = newName
    copied = type(doc)(
        projectName=doc.projectName,
        resource=doc.resource,
        dataFormat=doc.dataFormat,
        type=doc.type,
        desc=desc,
    )
    copied.save()
`;

// Duplicates the given documents, then reloads the project. Returns false on failure.
export const duplicateDocuments = async (docOids: string[]) => {
  const projectName = useProjectStore.getState().currProjectName;
  const { data } = await fetchPython({
    results: [],
    label: `duplicate ${docOids.length} documents`,
    code: [
      DUPLICATE_HELPER,
      ...docOids.map(oid => `duplicateDocumentByID('${oid}')`),
    ].join('\n'),
  });
  if (!data) {
    return false;
  }
  await fetchProjectDetails(projectName);
  return true;
};
