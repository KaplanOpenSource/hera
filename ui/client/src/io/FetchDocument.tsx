import { ProjectDocument } from "@shared/types";
import { fetchPython } from "./fetchPython";
import { FORBIDDEN_FIELDS } from "../shared/constants";

export const fetchDocument = async (docid: string) => {
  const { data } = await fetchPython({
    results: ['docData'],
    label: `get document ${docid}`,
    code: `
from hera.datalayer import All
docData = All.getDocumentByID('${docid}').asDict(with_id=True)
`,
  });
  return data?.docData;
}

const stringifyToPython = (obj: any) => {
  const jsonString = JSON.stringify(obj, (key, value) => {
    if (value === null) return '=-=None=-=';
    if (value === true) return '=-=True=-=';
    if (value === false) return '=-=False=-=';
    return value;
  });

  return jsonString
    .replace(/"=-=None=-="/g, 'None')
    .replace(/"=-=True=-="/g, 'True')
    .replace(/"=-=False=-="/g, 'False');
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
    label: 'update document',
    code: lines.join('\n'),
  });
  return data?.docData;
}
