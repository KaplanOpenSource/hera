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

const stringifyToPython = (obj: any) => {
  const jsonString = JSON.stringify(obj, (key, value) => {
    // Convert null to a unique placeholder for replacement
    return value === null ? '=-=None=-=' : value;
  });

  // Replace the placeholder 'None' without quotes
  return jsonString.replace(/"=-=None=-="/g, 'None');
}

export const updateDocument = async (newDoc: any, prevDoc: any) => {
  const docid = (prevDoc as ProjectDocument)._id.$oid;
  const lines = [`
from hera.datalayer import All
doc = All.getDocumentByID('${docid}')
`];
  for (const [field, prevVal] of Object.entries(prevDoc)) {
    if (!FORBIDDEN_FIELDS.includes(field) && JSON.stringify(prevVal) !== JSON.stringify(newDoc[field])) {
      lines.push(`doc.${field} = ${stringifyToPython(newDoc[field])}`)
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
