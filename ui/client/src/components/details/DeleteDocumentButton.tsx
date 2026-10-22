import { Delete } from '@mui/icons-material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { useConfirm } from '../../elements/useConfirm';
import { deleteDocuments } from '../../io/deleteDocuments';
import { ProjectDocument } from '../../shared/types';

// Deletes a single document from the document interface. The project's open tab
// closes on its own once the refreshed project no longer contains the document
// (see LayoutModel.closeMissingDocuments).
export const DeleteDocumentButton = ({
  document,
  displayName,
  isConfig,
}: {
  document: ProjectDocument,
  displayName: string,
  isConfig: boolean,
}) => {
  const { confirmOpen, ConfirmDialog } = useConfirm();

  return (
    <ButtonTooltip
      title={isConfig ? 'Project Config is deleted only with project' : 'Delete Document'}
      disabled={isConfig}
      onClick={async () => {
        if ((await confirmOpen({ title: `Are you sure you want to delete ${displayName}?` })).confirmed) {
          deleteDocuments([document._id.$oid]);
        }
      }}
    >
      <Delete />
      {ConfirmDialog}
    </ButtonTooltip>
  );
};
