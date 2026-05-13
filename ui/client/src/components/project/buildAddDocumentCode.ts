import { DocumentDesc } from "../../shared/types";
import { DocKind } from "./AddDocumentButton";
import { buildNotebookCommand } from "./buildNotebookCommand";

export const buildAddDocumentCode = ({
  kind,
  projectName,
  desc,
  toolkitNames,
  notebookResource,
  collection,
  resource,
}: {
  kind: DocKind,
  projectName: string,
  desc: DocumentDesc,
  toolkitNames: string[],
  notebookResource: string,
  collection: string,
  resource: string,
}) => {
  let addCommand: string;
  if (kind === DocKind.Notebook) {
    addCommand = buildNotebookCommand({ notebookResource, projectName, toolkitNames, desc });
  } else if (kind === DocKind.Agent) {
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
