import { Chip, TableCell, TableRow } from "@mui/material";
import { ProjectFolder } from "../../shared/types";
import { FilesFateEnum } from "./FilesFateEnum";
import { FolderStateNote } from "./FolderStateNote";

// One project in the delete confirm table: name, files folder, what happens to the files, and a folder note.
export const DeleteProjectsRow = ({
  project,
  fate,
}: {
  project: ProjectFolder,
  fate: FilesFateEnum,
}) => {
  return (
    <TableRow data-testid="delete-projects-row">
      <TableCell sx={{ overflowWrap: 'anywhere' }} data-testid="delete-projects-name">
        {project.name}
      </TableCell>
      <TableCell sx={{ overflowWrap: 'anywhere', color: 'text.secondary' }}>
        {project.filesDirectory}
      </TableCell>
      <TableCell data-testid="delete-projects-fate">
        {fate === FilesFateEnum.Removed && <Chip size="small" variant="outlined" label="Empty, removed" />}
        {fate === FilesFateEnum.Kept && <Chip size="small" color="success" variant="outlined" label="Kept" />}
        {fate === FilesFateEnum.Nothing && <Chip size="small" variant="outlined" label="None" />}
      </TableCell>
      <TableCell>
        <FolderStateNote state={project.folderState} />
      </TableCell>
    </TableRow>
  );
}
