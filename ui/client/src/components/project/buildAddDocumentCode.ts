import { DocumentDesc, Toolkit } from "../../shared/types";
import { buildNotebookCommand } from "./buildNotebookCommand";

export const buildAddDocumentCode = ({
  kind,
  projectName,
  desc,
  toolkits,
  notebookResource,
  collection,
  resource,
}: {
  kind: string,
  projectName: string,
  desc: DocumentDesc,
  toolkits: Toolkit[],
  notebookResource: string,
  collection: string,
  resource: string,
}) => {
  let addCommand: string;
  if (kind === 'Notebook') {
    addCommand = buildNotebookCommand({ notebookResource, projectName, toolkits, desc });
  } else if (kind === 'Agent') {
    addCommand = `
All.addDocument(
    '${projectName}',
    resource={"effects": {}},
    desc=${JSON.stringify(desc)},
    dataFormat=datatypes.JSON_DICT,
    type='ToolkitDataSource',
)`;
  } else {
    addCommand = `
${collection}().addDocument(
    '${projectName}',
    resource='${resource}',
    desc=${JSON.stringify(desc)},
)`;
  }

  return `
from hera.datalayer import All, datatypes, Measurements_Collection, Simulations_Collection, Cache_Collection
${addCommand}
docs = All.getDocumentsAsDict('${projectName}', with_id=True)
project = {"name": '${projectName}', "documents": docs['documents']}
`;
};
