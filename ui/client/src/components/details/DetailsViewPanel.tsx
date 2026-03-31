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
      <DetailsViewRepo
        repoPath={repoid}
      />
    )
  }

  const notebookName = idFromNotebookId(showItemId);
  if (notebookName) {
    return (
      <DetailsViewNotebook
        name={notebookName}
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

