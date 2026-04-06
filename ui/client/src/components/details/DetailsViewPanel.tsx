import { Box } from '@mui/material';
import { ProjectObj } from '../../objects/ProjectObj';
import { idFromDocId, idFromNotebookId, idFromRepoId } from '../../shared/idDocId';
import { DetailsViewDocId } from './DetailsViewDocId';
import { DetailsViewNotebook } from './DetailsViewNotebook';
import { DetailsViewRepo } from './DetailsViewRepo';

export const DetailsViewPanel = ({
  project,
  showItemId,
}: {
  project: ProjectObj,
  showItemId: string,
}) => {
  const docid = idFromDocId(showItemId);
  if (docid) {
    return (
      <DetailsViewDocId
        project={project}
        docid={docid}
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

  const notebookName = idFromNotebookId(showItemId);
  if (notebookName) {
    const filesDir = project?.configDocument?.data.desc.filesDirectory ?? '';
    return (
      <DetailsViewNotebook
        rootDir={filesDir}
        notebookName={notebookName}
      />
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

