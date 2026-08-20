import { useEffect, useRef, useState } from 'react';
import { DocumentObj } from '../../objects/ProjectObj';
import { ProjectDocument, WorkflowDesc } from '../../shared/types';
import { isWorkflowDoc } from '../../shared/workflow';
import { MutatorsListHandler } from '../../shared/workflowMutators/MutatorsListHandler';
import { DetailsViewDocumentContent } from './DetailsViewDocumentContent';

export const DetailsViewDocument = ({
  doc,
  setDoc,
}: {
  doc: DocumentObj,
  setDoc: (newDoc: DocumentObj) => Promise<void>,
}) => {
  const [shownDoc, setShownDoc] = useState<ProjectDocument>(JSON.parse(JSON.stringify(doc.data)));
  const lastLoadedRef = useRef(doc.data);

  // On reload, adopt the new server data only if the user has no unsaved edits;
  // if the document is being edited, keep their edits (don't reload it). This is
  // the load path, so it does NOT run the workflow mutators.
  useEffect(() => {
    setShownDoc(prev =>
      JSON.stringify(prev) === JSON.stringify(lastLoadedRef.current)
        ? JSON.parse(JSON.stringify(doc.data))
        : prev
    );
    lastLoadedRef.current = doc.data;
  }, [doc.data])

  // Every user edit goes through here; a workflow doc is run through its mutator pipeline before being stored.
  const changeShownDoc = (newDoc: ProjectDocument) => {
    if (isWorkflowDoc(newDoc)) {
      setShownDoc({ ...newDoc, desc: MutatorsListHandler.normalize(newDoc.desc as WorkflowDesc, doc.project.name) });
    } else {
      setShownDoc(newDoc);
    }
  };

  return (
    <DetailsViewDocumentContent
      doc={doc}
      setDoc={setDoc}
      shownDoc={shownDoc}
      setShownDoc={changeShownDoc}
    />
  );
}
