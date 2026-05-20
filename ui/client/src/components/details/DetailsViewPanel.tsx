import { Box } from '@mui/material';
import { ProjectObj } from '../../objects/ProjectObj';
import { CENTRAL_REPO_FOLDER_ID, idFromDocId, idFromRepoId } from '../../shared/idDocId';
import { DetailsViewDocId } from './DetailsViewDocId';
import { DetailsViewMergedRepo } from './DetailsViewMergedRepo';
import { DetailsViewRepo } from './DetailsViewRepo';

export const DetailsViewPanel = ({
  project,
  showItemId,
  previewHidden,
  setPreviewHidden,
}: {
  project: ProjectObj,
  showItemId: string,
  previewHidden: boolean,
  setPreviewHidden: (hidden: boolean) => void,
}) => {
  if (showItemId === CENTRAL_REPO_FOLDER_ID) {
    return (
      <Box sx={{ p: 2, height: '100%', overflow: 'auto' }}>
        <DetailsViewMergedRepo />
      </Box>
    )
  }

  const docid = idFromDocId(showItemId);
  if (docid) {
    return (
      <DetailsViewDocId
        project={project}
        docid={docid}
        previewHidden={previewHidden}
        setPreviewHidden={setPreviewHidden}
      />
    )
  }

  const repoid = idFromRepoId(showItemId);
  if (repoid) {
    return (
      <Box sx={{ p: 2, height: '100%', overflow: 'auto' }}>
        <DetailsViewRepo
          repoPath={repoid}
        />
      </Box>
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

