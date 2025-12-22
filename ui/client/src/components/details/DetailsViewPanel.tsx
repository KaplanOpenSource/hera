import { Paper } from '@mui/material';
import { useEffect, useState } from 'react';
import { fetchDocument, updateDocument } from '../../io/FetchDocument';
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
      if (docid) {
        const data = await fetchDocument(docid);
        if (data) (
          setDoc(data)
        )
      } else {
        setDoc(undefined);
      }
    })()
  }, [docid])

  const changeDocument = async (shownDoc: any) => {
    const data = await updateDocument(shownDoc, doc);
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
            setDoc={(newVal) => changeDocument(newVal)}
          />
        )
        : (
          <DetailsViewProject
            project={project}
          />)}
    </Paper>
  );
};
