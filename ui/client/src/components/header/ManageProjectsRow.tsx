import { TableCell, TableRow, Typography } from "@mui/material";
import { ProjectFolder } from "../../shared/types";
import { FolderCheckEnum } from "./FolderCheckEnum";
import { MarkProjectCheckbox } from "./MarkProjectCheckbox";

// One project in the manage table: a mark, its name, files folder, and an empty-folder note.
export const ManageProjectsRow = ({
  project,
  folderCheck,
  marked,
  onToggle,
}: {
  project: ProjectFolder,
  folderCheck: FolderCheckEnum,
  marked: boolean,
  onToggle: (shift: boolean) => void,
}) => {
  let folderText = 'No directory';
  if (project.filesDirectory) {
    folderText = project.filesDirectory;
  }

  return (
    <TableRow data-testid="manage-projects-row">
      <TableCell padding="checkbox">
        <MarkProjectCheckbox
          name={project.name}
          marked={marked}
          onToggle={onToggle}
        />
      </TableCell>
      <TableCell sx={{ overflowWrap: 'anywhere' }}>
        <span data-testid="manage-projects-name">{project.name}</span>
      </TableCell>
      <TableCell sx={{ overflowWrap: 'anywhere', color: 'text.secondary' }}>
        {folderText}
      </TableCell>
      <TableCell>
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
      </TableCell>
    </TableRow>
  );
}
