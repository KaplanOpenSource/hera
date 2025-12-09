import { Paper } from '@mui/material';
import { useEffect, useState } from 'react';
import { execPython } from '../../io/execPython';
import { idFromDocId } from '../../shared/idDocId';
import type { ProjectEntire } from '../../shared/types';
import { DetailsViewDocument } from './DetailsViewDocument';
import { DetailsViewProject } from './DetailsViewProject';

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
      console.log(docid);
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

  return (
    <Paper sx={{ p: 2, height: '100%' }}>
      {doc
        ? (
          <DetailsViewDocument
            doc={doc}
          />
        )
        : (
          <DetailsViewProject
            project={project}
          />)}
    </Paper>
  );
};
