import { DeleteForever, WarningAmber } from "@mui/icons-material";
import {
  Alert,
  Box,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  TextField,
  Typography
} from "@mui/material";
import { useState } from "react";
import { FolderStateEnum } from "../../shared/FolderStateEnum";
import { ProjectFolder } from "../../shared/types";
import { ExportProgress } from "./ExportProgress";
import { DeleteProjectsRow } from "./DeleteProjectsRow";
import { ExportStep } from "./ExportStep";
import { FilesFateEnum } from "./FilesFateEnum";
import { projectsText } from "./projectsText";

// What the delete does to a project's files folder, by its state.
const filesFate = (state: FolderStateEnum) => {
  const isEmpty = state === FolderStateEnum.Empty || state === FolderStateEnum.OnlyCaseConfig;
  const nothingOnDisk = state === FolderStateEnum.NoFolder || state === FolderStateEnum.Missing;
  if (isEmpty) {
    return FilesFateEnum.Removed;
  }
  if (nothingOnDisk) {
    return FilesFateEnum.Nothing;
  }
  return FilesFateEnum.Kept;
}

// Asks before deleting the marked projects, then stays open with the progress until the work ends.
export const DeleteProjectsConfirm = ({
  projects,
  running,
  step,
  onCancel,
  onConfirm,
}: {
  projects: ProjectFolder[],
  running: boolean,
  step: ExportStep | null,
  onCancel: () => void,
  onConfirm: () => void,
}) => {
  const [typed, setTyped] = useState('');

  // Rows with an empty folder come first, then kept, then nothing on disk.
  const rows: { project: ProjectFolder, fate: FilesFateEnum }[] = [];
  for (const fate of [FilesFateEnum.Removed, FilesFateEnum.Kept, FilesFateEnum.Nothing]) {
    for (const project of projects) {
      if (filesFate(project.folderState) === fate) {
        rows.push({ project, fate });
      }
    }
  }

  const count = projects.length;
  const countMatches = typed.trim() === String(count);
  const deleting = running && !step;

  // While the work runs the dialog cannot be closed.
  const close = () => {
    if (!running) {
      onCancel();
    }
  }

  let confirmLabel = `Export and delete ${projectsText(count)}`;
  if (running) {
    confirmLabel = 'Deleting...';
  }

  return (
    <Dialog open onClose={close} fullWidth maxWidth="md">
      <DialogTitle>
        <Stack direction="row" spacing={1} alignItems="center">
          <WarningAmber color="error" />
          <span>{`Delete ${projectsText(count)}?`}</span>
        </Stack>
      </DialogTitle>
      <DialogContent dividers>
        <TableContainer sx={{ maxHeight: 320 }}>
          <Table size="small" stickyHeader sx={{ tableLayout: 'fixed' }}>
            <TableHead>
              <TableRow>
                <TableCell sx={{ width: '22%' }}>Name</TableCell>
                <TableCell sx={{ width: '33%' }}>Files folder</TableCell>
                <TableCell sx={{ width: '20%' }}>Folder</TableCell>
                <TableCell sx={{ width: '25%' }}>Note</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {rows.map((row) => (
                <DeleteProjectsRow
                  key={row.project.name}
                  project={row.project}
                  fate={row.fate}
                />
              ))}
            </TableBody>
          </Table>
        </TableContainer>
        <Alert severity="info" variant="filled" sx={{ my: 2 }}>
          projects.zip is saved first. If any export fails, nothing is deleted.
          Files are never deleted. Only an empty folder is removed.
        </Alert>
        <TextField
          autoFocus
          fullWidth
          size="small"
          label={`Type ${count} to confirm`}
          value={typed}
          onChange={(e) => setTyped(e.target.value)}
          disabled={running}
          slotProps={{ htmlInput: { 'data-testid': 'delete-projects-count' } }}
        />
        <Box sx={{ mt: 2 }}>
          {step && <ExportProgress step={step} />}
          {deleting && (
            <Typography variant="body2" color="text.secondary" data-testid="delete-projects-deleting">
              Deleting the projects...
            </Typography>
          )}
        </Box>
      </DialogContent>
      <DialogActions sx={{ px: 3, py: 1.5 }}>
        <Button onClick={onCancel} disabled={running}>
          Cancel
        </Button>
        <Button
          variant="contained"
          color="error"
          startIcon={<DeleteForever />}
          disabled={!countMatches || running}
          onClick={onConfirm}
          data-testid="delete-projects-confirm"
        >
          {confirmLabel}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
