import { Paper } from '@mui/material';
import { useEffect, useState } from 'react';
import { fetchDocument, updateDocument } from '../../io/FetchDocument';
import { DocumentObj, ProjectObj } from '../../objects/ProjectObj';
import { idFromDocId } from '../../shared/idDocId';
import { DetailsViewDocument } from './DetailsViewDocument';

export const DetailsViewPanel = ({
  project,
  showItemId,
}: {
  project: ProjectObj,
  showItemId: string,
}) => {
  const [doc, setDoc] = useState<any>(undefined);

  const docid = idFromDocId(showItemId);
  if (docid) {
    return (
      <DetailsViewDocId
        project={project}
        docid={docid}
      />
    )
  }

  if (project?.configDocument?.docid) {
    return (
      <DetailsViewDocId
        project={project}
        docid={project?.configDocument?.docid}
      />
    )
  } else {
    return null;
  }
};


export const DetailsViewDocId = ({
  project,
  docid,
}: {
  project: ProjectObj,
  docid: string,
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
    })()
  }, [docid, project?.name])

  const changeDocument = async (shownDoc: any) => {
    const data = await updateDocument(shownDoc, doc);
    if (data) {
      setDoc(data)
    }
  }

  return (
    <>
      {doc
        ? (
          <DetailsViewDocument
            doc={new DocumentObj(doc, project)}
            setDoc={(newDoc) => changeDocument(newDoc.data)}
          />
        )
        : null}
    </>
  );
};
