import { useMemo } from 'react';
import { updateDocument } from '../../io/FetchDocument';
import { fetchProjectDetails } from '../../io/FetchProjects';
import { DocumentObj } from '../../objects/ProjectObj';
import { useProjectStore } from '../../stores/useProjectStore';
import { ProjectDocument, WorkflowDesc } from '../../shared/types';
import { isWorkflowDoc } from '../../shared/workflow';
import { DocumentFieldsMutator } from '../../shared/workflowMutators/DocumentFieldsMutator';

// The document as the user sees it: their unsaved edits if there are any, else a
// copy of what was loaded. Shared by every view of the same document.
export const useShownDoc = (doc: DocumentObj) => {
  const editedDoc = useProjectStore(state => state.editedDocs[doc.docid]);
  const setEditedDoc = useProjectStore(state => state.setEditedDoc);

  // A copy, so nothing editing the shown document can touch the stored one.
  const loadedDoc = useMemo(() => JSON.parse(JSON.stringify(doc.data)), [doc.data]);

  const shownDoc: ProjectDocument = editedDoc ?? loadedDoc;

  // Every user edit goes through here; a workflow doc is run through its mutator pipeline before being stored.
  const setShownDoc = (newDoc: ProjectDocument) => {
    if (isWorkflowDoc(newDoc)) {
      setEditedDoc(doc.docid, { ...newDoc, desc: DocumentFieldsMutator.mutate(newDoc.desc as WorkflowDesc) });
    } else {
      setEditedDoc(doc.docid, newDoc);
    }
  };

  // Stores the edits, then reloads the project so every view shows the saved value.
  const saveShownDoc = async () => {
    if (!editedDoc) {
      return;
    }
    const data = await updateDocument(editedDoc, doc.data);
    if (data) {
      await fetchProjectDetails(doc.project.name, true);
      // Dropped only now, so the view never falls back to the pre-save document.
      setEditedDoc(doc.docid, null);
    }
  };

  return { shownDoc, setShownDoc, saveShownDoc };
};
