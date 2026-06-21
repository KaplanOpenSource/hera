import { Box } from '@mui/material';
import { ProjectObj } from '../../objects/ProjectObj';
import { classifyItemId, idFromDocId, idFromRepoId, ItemKind, toolkitNameFromSplitId } from '../../shared/idDocId';
import { DetailsViewDocId } from './DetailsViewDocId';
import { DetailsViewMergedRepo } from './DetailsViewMergedRepo';
import { DetailsViewRepo } from './DetailsViewRepo';
import { DetailsViewToolkit } from './toolkit/DetailsViewToolkit';

export { detailsTabName } from '../../shared/tabKind';

export const DetailsViewPanel = ({
  project,
  showItemId,
}: {
  project: ProjectObj,
  showItemId: string,
}) => {
  const kind = classifyItemId(showItemId);

  if (kind === ItemKind.CentralRepo) {
    return (
      <Box sx={{ p: 2, height: '100%', overflow: 'auto' }}>
        <DetailsViewMergedRepo />
      </Box>
    )
  }

  if (kind === ItemKind.Split) {
    const toolkitName = toolkitNameFromSplitId(showItemId, project.documents);
    if (toolkitName) {
      return (
        <Box sx={{ height: '100%', overflow: 'auto' }}>
          <DetailsViewToolkit toolkitName={toolkitName} />
        </Box>
      )
    }
  }

  if (kind === ItemKind.Document) {
    const docid = idFromDocId(showItemId);
    if (docid) {
      return (
        <DetailsViewDocId
          project={project}
          docid={docid}
        />
      )
    }
  }

  if (kind === ItemKind.Repo) {
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
