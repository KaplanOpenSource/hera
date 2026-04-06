import { Box } from '@mui/material';
import { useEffect, useState } from 'react';
import { fetchDocument, updateDocument } from '../../io/FetchDocument';
import { DocumentObj, ProjectObj } from '../../objects/ProjectObj';
import { DetailsViewDocument } from './DetailsViewDocument';
import { DetailsViewNotebook } from './DetailsViewNotebook';

export const DetailsViewDocId = ({
  project, docid,
}: {
  project: ProjectObj;
  docid: string;
}) => {
  const [doc, setDoc] = useState<any>(undefined);

  useEffect(() => {
    (async () => {
      if (docid) {
        const data = await fetchDocument(docid);
        if (data) {
          setDoc(data);
          return;
        }
      }
      setDoc(undefined);
    })();
  }, [docid, project]);

  const changeDocument = async (shownDoc: any) => {
    const data = await updateDocument(shownDoc, doc);
    if (data) {
      setDoc(data);
    }
  };

  const docObj = doc ? new DocumentObj(doc, project) : null;

  return (
    <>
      {docObj
        ? docObj.isNotebook
          ? (
            <DetailsViewNotebook
              rootDir={project.configDocument?.data.desc.filesDirectory ?? ''}
              notebookName={docObj.data.desc.datasourceName ?? ''}
            />
          )
          : (
            <Box sx={{ p: 2, height: '100%', overflow: 'auto' }}>
              <DetailsViewDocument
                doc={docObj}
                setDoc={(newDoc) => changeDocument(newDoc.data)}
              />
            </Box>
          )
        : null}
    </>
  );
};
