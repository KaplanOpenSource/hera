import { ProjectObj } from '../../objects/ProjectObj';
import { idFromDocId } from '../../shared/idDocId';
import { DetailsViewDocId } from './DetailsViewDocId';

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

