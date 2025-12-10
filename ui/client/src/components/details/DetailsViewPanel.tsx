import { Paper } from '@mui/material';
import { useEffect, useState } from 'react';
import { execPython } from '../../io/execPython';
import { idFromDocId } from '../../shared/idDocId';
import type { ProjectDocument, ProjectEntire } from '../../shared/types';
import { DetailsViewDocument } from './DetailsViewDocument';
import { DetailsViewProject } from './DetailsViewProject';

const FORBIDDEN_FIELDS = ['_id', '_cls', 'projectName'];

export const DetailsViewPanel = ({
  project,
  selectedItemsIds,
}: {
  project: ProjectEntire,
  selectedItemsIds: string[],
}) => {
  const docid = idFromDocId(selectedItemsIds[0]);
  const [doc, setDoc] = useState<any>(undefined);

  useEffect(() => {
    (async () => {
      if (docid) {
        const { data } = await execPython(`
import json
from hera.datalayer import All
docs = All.getDocumentByID('${docid}')
result = docs.asDict(with_id=True)
`);
        if (data) (
          setDoc(data)
        )
      } else {
        setDoc(undefined);
      }
    })()
  }, [docid])


  const updateDocument = async (shownDoc: any) => {
    const docid = (doc as ProjectDocument)._id.$oid;
    const lines = [`
from hera.datalayer import All
doc = All.getDocumentByID('${docid}')
`];
    for (const [field, prevVal] of Object.entries(doc)) {
      if (!FORBIDDEN_FIELDS.includes(field) && JSON.stringify(prevVal) !== JSON.stringify(shownDoc[field])) {
        lines.push(`doc.${field} = ${JSON.stringify(shownDoc[field])}`)
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
      setDoc(data)
    }
  }

  return (
    <Paper sx={{ p: 2, height: '100%' }}>
      {doc
        ? (
          <DetailsViewDocument
            doc={doc}
            setDoc={(newVal) => updateDocument(newVal)}
          />
        )
        : (
          <DetailsViewProject
            project={project}
          />)}
    </Paper>
  );
};
