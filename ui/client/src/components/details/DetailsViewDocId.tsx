import { useEffect, useState } from 'react';
import { fetchDocument, updateDocument } from '../../io/FetchDocument';
import { DocumentObj, ProjectObj } from '../../objects/ProjectObj';
import { DetailsViewDocument } from './DetailsViewDocument';

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
  }, [docid, project?.name]);

  const changeDocument = async (shownDoc: any) => {
    const data = await updateDocument(shownDoc, doc);
    if (data) {
      setDoc(data);
    }
  };

  return (
    <>
      {doc
        ? (
          <DetailsViewDocument
            doc={new DocumentObj(doc, project)}
            setDoc={(newDoc) => changeDocument(newDoc.data)} />
        )
        : null}
    </>
  );
};
