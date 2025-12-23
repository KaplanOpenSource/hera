import { Paper } from '@mui/material';
import { useEffect, useState } from 'react';
import { fetchDocument, updateDocument } from '../../io/FetchDocument';
import { DocumentObj, ProjectObj } from '../../objects/ProjectObj';
import { idFromDocId } from '../../shared/idDocId';
import { DetailsViewDocument } from './DetailsViewDocument';

export const DetailsViewPanel = ({
  project,
  selectedItemsIds,
}: {
  project: ProjectObj,
  selectedItemsIds: string[],
}) => {
  const [doc, setDoc] = useState<any>(undefined);

  useEffect(() => {
    (async () => {
      const docid = idFromDocId(selectedItemsIds[0]);
      if (docid) {
        const data = await fetchDocument(docid);
        if (data) {
          setDoc(data);
          return;
        }
      }
      const confid = project?.configDocument?.docid;
      if (confid) {
        const data = await fetchDocument(confid);
        if (data) {
          setDoc(data);
          return;
        }
      }
      setDoc(undefined);
    })()
  }, [selectedItemsIds[0]])

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
            doc={new DocumentObj(doc, project)}
            setDoc={(newDoc) => changeDocument(newDoc.data)}
          />
        )
        : null}
    </Paper>
  );
};
