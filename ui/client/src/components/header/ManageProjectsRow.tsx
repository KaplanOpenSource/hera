import { ListItem, ListItemText, Stack, Typography } from "@mui/material";
import { ProjectFolder } from "../../shared/types";
import { FolderCheckEnum } from "./FolderCheckEnum";

// One project in the manage list: its name, files folder, and an empty-folder note.
export const ManageProjectsRow = ({
  project,
  folderCheck,
}: {
  project: ProjectFolder,
  folderCheck: FolderCheckEnum,
}) => {
  return (
    <ListItem divider data-testid="manage-projects-row">
      <ListItemText
        primary={
          <Stack direction="row" alignItems="center" spacing={1}>
            <span data-testid="manage-projects-name">{project.name}</span>
            {folderCheck === FolderCheckEnum.Checking && (
              <Typography variant="caption" color="text.secondary" data-testid="manage-projects-checking">
                checking...
              </Typography>
            )}
            {folderCheck === FolderCheckEnum.Missing && (
              <Typography variant="caption" color="warning.main" data-testid="manage-projects-missing">
                folder missing
              </Typography>
            )}
            {folderCheck === FolderCheckEnum.Empty && (
              <Typography variant="caption" color="warning.main" data-testid="manage-projects-empty">
                empty folder
              </Typography>
            )}
          </Stack>
        }
        secondary={project.filesDirectory || 'No directory'}
      />
    </ListItem>
  );
}
