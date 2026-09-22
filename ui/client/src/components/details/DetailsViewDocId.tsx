import { Box } from '@mui/material';
import { ProjectObj } from '../../objects/ProjectObj';
import { DetailsViewDocumentContent } from './DetailsViewDocumentContent';
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
              <DetailsViewDocumentContent doc={docObj} />
            </Box>
          )
        : null}
    </>
  );
};
