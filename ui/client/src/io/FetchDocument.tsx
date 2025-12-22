import { ProjectDocument } from "@shared/types";
import { execPython } from "./execPython";

const FORBIDDEN_FIELDS = ['_id', '_cls', 'projectName'];

export const fetchDocument = async (docid: string) => {
  const { data } = await execPython(`
import json
from hera.datalayer import All
docs = All.getDocumentByID('${docid}')
result = docs.asDict(with_id=True)
`);
  if (data) {
    return data;
  }
  return undefined;
}

export const updateDocument = async (newDoc: any, prevDoc: any) => {
  const docid = (prevDoc as ProjectDocument)._id.$oid;
  const lines = [`
from hera.datalayer import All
doc = All.getDocumentByID('${docid}')
`];
  for (const [field, prevVal] of Object.entries(prevDoc)) {
    if (!FORBIDDEN_FIELDS.includes(field) && JSON.stringify(prevVal) !== JSON.stringify(newDoc[field])) {
      lines.push(`doc.${field} = ${JSON.stringify(newDoc[field])}`)
    }
  }
  lines.push(`
doc.save()
docs = All.getDocumentByID('${docid}')
result = docs.asDict(with_id=True)
`)
  const code = lines.join('\n');
  const { data } = await execPython(code);
  if (data) {
    return data;
  }
  return undefined;
}
