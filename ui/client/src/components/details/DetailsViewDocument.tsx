import { useEffect, useState } from 'react';
import { DocumentObj } from '../../objects/ProjectObj';
import { ProjectDocument } from '../../shared/types';
import { DetailsViewDocumentContent } from './DetailsViewDocumentContent';

export const DetailsViewDocument = ({
  doc,
  setDoc,
}: {
  doc: DocumentObj,
  setDoc: (newDoc: DocumentObj) => void,
}) => {
  const [shownDoc, setShownDoc] = useState<ProjectDocument>(JSON.parse(JSON.stringify(doc.data)));

  useEffect(() => {
    setShownDoc(JSON.parse(JSON.stringify(doc.data)));
  }, [doc.data])

  return (
    <DetailsViewDocumentContent
      doc={doc}
      setDoc={setDoc}
      shownDoc={shownDoc}
      setShownDoc={setShownDoc}
    />
  );
}
