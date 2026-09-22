import { Box } from '@mui/material';
import { updateDocument } from '../../io/FetchDocument';
import { fetchProjectDetails } from '../../io/FetchProjects';
import { ProjectObj } from '../../objects/ProjectObj';
import { useProjectStore } from '../../stores/useProjectStore';
import { DetailsViewDocument } from './DetailsViewDocument';
import { DetailsViewNotebook } from './DetailsViewNotebook';

export const DetailsViewDocId = ({
  project,
  docid,
}: {
  project: ProjectObj;
  docid: string;
}) => {
  // The document data comes from the project store (loaded centrally and auto-reloaded),
  // so there is no per-tab fetch — open tabs stay in sync with the one store.
  const docObj = project.allDocuments.find(d => d.docid === docid) ?? null;
  const setEditedDoc = useProjectStore(state => state.setEditedDoc);

  const changeDocument = async (shownDoc: any) => {
    if (!docObj) return;
    const data = await updateDocument(shownDoc, docObj.data);
    if (data) {
      // Pull the saved value back into the store so all views update.
      await fetchProjectDetails(project.name, true);
      // Dropped only now, so the view never falls back to the pre-save document.
      setEditedDoc(docid, null);
    }
  };

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
            <Box sx={{ p: 2, height: '100%', overflow: 'auto', display: 'flex', flexDirection: 'column', minHeight: 0 }}>
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
