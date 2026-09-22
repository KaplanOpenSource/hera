import { useMemo } from 'react';
import { DocumentObj } from '../../objects/ProjectObj';
import { useProjectStore } from '../../stores/useProjectStore';
import { ProjectDocument, WorkflowDesc } from '../../shared/types';
import { isWorkflowDoc } from '../../shared/workflow';
import { DocumentFieldsMutator } from '../../shared/workflowMutators/DocumentFieldsMutator';
import { DetailsViewDocumentContent } from './DetailsViewDocumentContent';

export const DetailsViewDocument = ({
  doc,
  setDoc,
}: {
  doc: DocumentObj,
  setDoc: (newDoc: DocumentObj) => Promise<void>,
}) => {
  const editedDoc = useProjectStore(state => state.editedDocs[doc.docid]);
  const setEditedDoc = useProjectStore(state => state.setEditedDoc);

  // A copy, so nothing editing the shown document can touch the stored one.
  const loadedDoc = useMemo(() => JSON.parse(JSON.stringify(doc.data)), [doc.data]);

  // While the document has unsaved edits, show those; a reload leaves them alone.
  const shownDoc = editedDoc ?? loadedDoc;

  // Every user edit goes through here; a workflow doc is run through its mutator pipeline before being stored.
  const changeShownDoc = (newDoc: ProjectDocument) => {
    if (isWorkflowDoc(newDoc)) {
      setEditedDoc(doc.docid, { ...newDoc, desc: DocumentFieldsMutator.mutate(newDoc.desc as WorkflowDesc) });
    } else {
      setEditedDoc(doc.docid, newDoc);
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
