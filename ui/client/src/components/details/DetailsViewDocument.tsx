import { useEffect, useRef, useState } from 'react';
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
  const lastLoadedRef = useRef(doc.data);

  // On reload, adopt the new server data only if the user has no unsaved edits;
  // if the document is being edited, keep their edits (don't reload it).
  useEffect(() => {
    setShownDoc(prev =>
      JSON.stringify(prev) === JSON.stringify(lastLoadedRef.current)
        ? JSON.parse(JSON.stringify(doc.data))
        : prev
    );
    lastLoadedRef.current = doc.data;
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
