import { Close, ContentCopy, Delete } from '@mui/icons-material';
import { Divider, ListItemIcon, ListItemText, Menu, MenuItem } from '@mui/material';
import { useConfirm } from '../../elements/useConfirm';
import { deleteDocuments } from '../../io/deleteDocuments';
import { duplicateDocuments } from '../../io/duplicateDocuments';

export type MenuPosition = {
  x: number,
  y: number,
};

export type MenuDocument = {
  oid: string,
  name: string,
};

// Right-click menu of the documents tree: acts on all the selected documents.
export const TreeContextMenu = ({
  position,
  docs,
  onClose,
  onDeleted,
}: {
  position: MenuPosition | null,
  docs: MenuDocument[],
  onClose: () => void,
  onDeleted: () => void,
}) => {
  const { confirmOpen, ConfirmDialog } = useConfirm();
  const docOids = docs.map(d => d.oid);
  const plural = docOids.length === 1 ? '' : 's';

  // One document is named in the question, many are counted.
  let deleteTitle = `Delete ${docOids.length} documents?`;
  if (docs.length === 1) {
    deleteTitle = `Delete ${docs[0].name}?`;
  }

  const handleDelete = async () => {
    onClose();
    const { confirmed } = await confirmOpen({
      title: deleteTitle,
    });
    if (!confirmed) return;
    const deleted = await deleteDocuments(docOids);
    if (deleted) {
      onDeleted();
    }
  };

  const handleDuplicate = async () => {
    onClose();
    await duplicateDocuments(docOids);
  };

  let anchorPosition = undefined;
  if (position) {
    anchorPosition = { top: position.y, left: position.x };
  }

  return (
    <>
      <Menu
        open={Boolean(position)}
        anchorReference="anchorPosition"
        anchorPosition={anchorPosition}
        onClose={onClose}
      >
        <MenuItem disabled={docOids.length === 0} onClick={handleDuplicate}>
          <ListItemIcon>
            <ContentCopy fontSize="small" />
          </ListItemIcon>
          <ListItemText>
            Duplicate {docOids.length} document{plural}
          </ListItemText>
        </MenuItem>
        <MenuItem disabled={docOids.length === 0} onClick={handleDelete}>
          <ListItemIcon>
            <Delete fontSize="small" color="error" />
          </ListItemIcon>
          <ListItemText>
            Delete {docOids.length} document{plural}
          </ListItemText>
        </MenuItem>
        <Divider />
        <MenuItem onClick={onClose}>
          <ListItemIcon>
            <Close fontSize="small" />
          </ListItemIcon>
          <ListItemText>
            Cancel
          </ListItemText>
        </MenuItem>
      </Menu>
      {ConfirmDialog}
    </>
  );
};
