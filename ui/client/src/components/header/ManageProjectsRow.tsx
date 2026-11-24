import { TableCell, TableRow } from "@mui/material";
import { ProjectFolder } from "../../shared/types";
import { FolderStateNote } from "./FolderStateNote";
import { MarkProjectCheckbox } from "./MarkProjectCheckbox";
import { SwitchToProjectButton } from "./SwitchToProjectButton";

// One project in the manage table: a mark, its name, files folder, document count, a folder note, and a switch button.
export const ManageProjectsRow = ({
  project,
  marked,
  onToggle,
  onSwitched,
}: {
  project: ProjectFolder,
  marked: boolean,
  onToggle: (shift: boolean) => void,
  onSwitched: () => void,
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
      <TableCell data-testid="manage-projects-documents">
        {project.documentCount}
      </TableCell>
      <TableCell>
        <FolderStateNote state={project.folderState} />
      </TableCell>
      <TableCell>
        <SwitchToProjectButton
          name={project.name}
          onSwitched={onSwitched}
        />
      </TableCell>
    </TableRow>
  );
}
