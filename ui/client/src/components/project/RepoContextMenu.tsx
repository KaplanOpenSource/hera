import { Close, DriveFolderUpload } from '@mui/icons-material';
import { Divider, ListItemIcon, ListItemText, Menu, MenuItem } from '@mui/material';
import { MenuPosition } from './TreeContextMenu';

// Right-click menu of a repository row: the same action as a double click.
export const RepoContextMenu = ({
  position,
  repoName,
  onClose,
  onLoad,
}: {
  position: MenuPosition | null,
  repoName: string | undefined,
  onClose: () => void,
  onLoad: (repoName: string) => void,
}) => {
  const handleLoad = async () => {
    onClose();
    if (repoName) {
      await onLoad(repoName);
    }
  };

  let anchorPosition = undefined;
  if (position) {
    anchorPosition = { top: position.y, left: position.x };
  }

  return (
    <Menu
      open={Boolean(position)}
      anchorReference="anchorPosition"
      anchorPosition={anchorPosition}
      onClose={onClose}
    >
      <MenuItem disabled={!repoName} onClick={handleLoad}>
        <ListItemIcon>
          <DriveFolderUpload fontSize="small" />
        </ListItemIcon>
        <ListItemText>
          Load "{repoName}" into project
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
  );
};
