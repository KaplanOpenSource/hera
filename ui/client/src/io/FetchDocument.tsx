import { ProjectDocument } from "@shared/types";
import { fetchPython } from "./fetchPython";
import { FORBIDDEN_FIELDS } from "../shared/constants";

export const fetchDocument = async (docid: string) => {
  const { data } = await fetchPython({
    results: ['docData'],
    code: `
from hera.datalayer import All
docData = All.getDocumentByID('${docid}').asDict(with_id=True)
`,
  });
  return data?.docData;
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
docData = All.getDocumentByID('${docid}').asDict(with_id=True)
`)
  const { data } = await fetchPython({
    results: ['docData'],
    code: lines.join('\n'),
  });
  return data?.docData;
}
