import { Box } from '@mui/material';
import { updateDocument } from '../../io/FetchDocument';
import { fetchProjectDetails } from '../../io/FetchProjects';
import { ProjectObj } from '../../objects/ProjectObj';
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

  const changeDocument = async (shownDoc: any) => {
    if (!docObj) return;
    const data = await updateDocument(shownDoc, docObj.data);
    if (data) {
      // Pull the saved value back into the store so all views update.
      fetchProjectDetails(project.name, true);
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
