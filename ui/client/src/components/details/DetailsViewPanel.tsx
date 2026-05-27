import { Box } from '@mui/material';
import { ProjectObj } from '../../objects/ProjectObj';
import { CENTRAL_REPO_FOLDER_ID, idFromDocId, idFromRepoId, idFromToolkitSplitId } from '../../shared/idDocId';
import { DetailsViewDocId } from './DetailsViewDocId';
import { DetailsViewMergedRepo } from './DetailsViewMergedRepo';
import { DetailsViewRepo } from './DetailsViewRepo';
import { DetailsViewToolkit } from './DetailsViewToolkit';

export const detailsTabName = (showItemId: string, project: ProjectObj): string => {
  if (showItemId === CENTRAL_REPO_FOLDER_ID) return 'Repositories';
  const toolkitName = idFromToolkitSplitId(showItemId);
  if (toolkitName) return toolkitName;
  const docid = idFromDocId(showItemId);
  if (docid) {
    const doc = project.allDocuments.find(d => d.docid === docid);
    if (doc) return doc.isConfig ? project.name + ' config' : doc.name;
  }
  const repoPath = idFromRepoId(showItemId);
  if (repoPath) return repoPath;
  return project.name + ' config';
};

export const DetailsViewPanel = ({
  project,
  showItemId,
}: {
  project: ProjectObj,
  showItemId: string,
}) => {
  if (showItemId === CENTRAL_REPO_FOLDER_ID) {
    return (
      <Box sx={{ p: 2, height: '100%', overflow: 'auto' }}>
        <DetailsViewMergedRepo />
      </Box>
    )
  }

  const toolkitName = idFromToolkitSplitId(showItemId);
  if (toolkitName) {
    return (
      <Box sx={{ height: '100%', overflow: 'auto' }}>
        <DetailsViewToolkit toolkitName={toolkitName} />
      </Box>
    )
  }

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

