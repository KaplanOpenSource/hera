import { Delete, Description } from '@mui/icons-material';
import { Stack, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { useConfirm } from '../../elements/useConfirm';
import { idNotebookId } from '../../shared/idDocId';

export const NotebookListItem = ({
  name,
  selectedItemId,
  onDelete,
  onRefresh,
}: {
  name: string,
  selectedItemId: string | undefined,
  onDelete: (name: string) => void,
  onRefresh: () => void,
}) => {
  const { confirmOpen, ConfirmDialog } = useConfirm();

  return (
    <TreeItem
      itemId={idNotebookId(name)}
      label={(
        <Stack direction="row" alignItems="center" onClick={() => {
          if (selectedItemId === idNotebookId(name)) onRefresh();
        }}>
          <Typography marginRight={1}>{name}</Typography>
          <ButtonTooltip title="Delete" onClick={async () => {
            if ((await confirmOpen({ title: `Are you sure you want to delete ${name}?` })).confirmed) {
              onDelete(name);
            }
          }}>
            <Delete fontSize="small" />
          </ButtonTooltip>
          {ConfirmDialog}
        </Stack>
      )}
      slots={{ icon: () => <Description fontSize="small" /> }}
    />
  );
};
