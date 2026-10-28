import { DeleteSweep } from '@mui/icons-material';
import { Stack } from '@mui/material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { useConfirm } from '../../elements/useConfirm';
import { deleteDocuments } from '../../io/deleteDocuments';
import { idFromDocId } from '../../shared/idDocId';
import { useProjectStore } from '../../stores/useProjectStore';

// Bulk-deletes the documents among the currently selected tree nodes.
export const DeleteSelectedButton = ({
  selectedIds = [],
  onDeleted,
}: {
  selectedIds?: string[],
  onDeleted?: () => void,
}) => {
  const { confirmOpen, ConfirmDialog } = useConfirm();

  // Of the selected tree nodes, keep the documents (excluding the project config doc).
  const configDocId = useProjectStore.getState().getProject()?.configDocument?.docid;
  const selectedDocOids = selectedIds
    .map(idFromDocId)
    .filter((oid): oid is string => !!oid && oid !== configDocId);

  const deleteSelected = async () => {
    const { confirmed } = await confirmOpen({
      title: `Delete ${selectedDocOids.length} selected document${selectedDocOids.length === 1 ? '' : 's'}?`,
    });
    if (!confirmed) return;

    const deleted = await deleteDocuments(selectedDocOids);
    if (!deleted) return;
    onDeleted?.();
  };

  return (
    <>
      <ButtonTooltip
        button
        color="error"
        title={`Delete selected (${selectedDocOids.length})`}
        disabled={selectedDocOids.length === 0}
        onClick={deleteSelected}
      >
        <Stack direction={'row'} alignItems={'center'} spacing={0.5}>
          <DeleteSweep fontSize="small" sx={{ transform: 'translateY(-1px)' }} />
        </Stack>
      </ButtonTooltip>
      {ConfirmDialog}
    </>
  );
};
