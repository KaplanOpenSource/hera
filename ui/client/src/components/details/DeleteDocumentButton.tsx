import { Delete } from '@mui/icons-material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { useConfirm } from '../../elements/useConfirm';
import { fetchProjectDetails } from '../../io/FetchProjects';
import { fetchPython } from '../../io/fetchPython';
import { ProjectDocument } from '../../shared/types';

// Deletes a single document from the document interface. The project's open tab
// closes on its own once the refreshed project no longer contains the document
// (see LayoutModel.closeMissingDocuments).
export const DeleteDocumentButton = ({
  document,
  projectName,
  displayName,
  isConfig,
}: {
  document: ProjectDocument,
  projectName: string,
  displayName: string,
  isConfig: boolean,
}) => {
  const { confirmOpen, ConfirmDialog } = useConfirm();

  const deleteDocument = async () => {
    const { data } = await fetchPython({
      results: [],
      label: `delete document ${displayName}`,
      code: `
from hera.datalayer import All
All.deleteDocumentByID('${document?._id.$oid}')
`,
    });
    if (!data) {
      return;
    }
    await fetchProjectDetails(projectName);
  };

  return (
    <ButtonTooltip
      title={isConfig ? 'Project Config is deleted only with project' : 'Delete Document'}
      disabled={isConfig}
      onClick={async () => {
        if ((await confirmOpen({ title: `Are you sure you want to delete ${displayName}?` })).confirmed) {
          deleteDocument();
        }
      }}
    >
      <Delete />
      {ConfirmDialog}
    </ButtonTooltip>
  );
};
