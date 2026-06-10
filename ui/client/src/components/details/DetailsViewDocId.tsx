import { Box } from '@mui/material';
import { useEffect, useState } from 'react';
import { fetchDocument, updateDocument } from '../../io/FetchDocument';
import { DocumentObj, ProjectObj } from '../../objects/ProjectObj';
import { DetailsViewDocument } from './DetailsViewDocument';
import { DetailsViewNotebook } from './DetailsViewNotebook';

export const DetailsViewDocId = ({
  project,
  docid,
}: {
  project: ProjectObj;
  docid: string;
}) => {
  const [doc, setDoc] = useState<any>(undefined);

  useEffect(() => {
    // Periodic reloads are silent (no notification); the initial open is not.
    const load = async (silent: boolean) => {
      // Skip docs that no longer exist (e.g. just deleted) to avoid a failing fetch.
      if (docid && project.documentIds.has(docid)) {
        const data = await fetchDocument(docid, silent);
        if (data) {
          setDoc(data);
          return;
        }
      }
      setDoc(undefined);
    };
    load(false);
    const interval = setInterval(() => load(true), 5000);
    return () => clearInterval(interval);
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
              resource={docObj.data.resource as string}
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
