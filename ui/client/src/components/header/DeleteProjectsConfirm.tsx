import {
  Box,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  TextField,
  Typography
} from "@mui/material";
import { useState } from "react";
import { ProjectFolder } from "../../shared/types";
import { ExportProgress } from "./ExportProgress";
import { ExportStep } from "./ExportStep";
import { FolderCheckEnum } from "./FolderCheckEnum";
import { projectsText } from "./projectsText";

// Asks before deleting the marked projects, then stays open with the progress until the work ends.
export const DeleteProjectsConfirm = ({
  projects,
  folderChecks,
  running,
  step,
  onCancel,
  onConfirm,
}: {
  projects: ProjectFolder[],
  folderChecks: { [name: string]: FolderCheckEnum },
  running: boolean,
  step: ExportStep | null,
  onCancel: () => void,
  onConfirm: () => void,
}) => {
  const [typed, setTyped] = useState('');

  const emptyFolder: ProjectFolder[] = [];
  const withFiles: ProjectFolder[] = [];
  const noFolder: ProjectFolder[] = [];
  for (const project of projects) {
    const check = folderChecks[project.name] ?? FolderCheckEnum.Checking;
    if (check === FolderCheckEnum.Empty) {
      emptyFolder.push(project);
    } else if (check === FolderCheckEnum.NoFolder || check === FolderCheckEnum.Missing) {
      noFolder.push(project);
    } else {
      withFiles.push(project);
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
    <Dialog open onClose={close} fullWidth maxWidth="sm">
      <DialogTitle>{`Delete ${projectsText(count)}?`}</DialogTitle>
      <DialogContent>
        {emptyFolder.length > 0 && (
          <Box sx={{ mb: 2 }} data-testid="delete-projects-empty">
            <Typography variant="subtitle2">Empty folder: their files are deleted too.</Typography>
            {emptyFolder.map((project) => (
              <Typography key={project.name} variant="body2">{project.name}</Typography>
            ))}
          </Box>
        )}
        {withFiles.length > 0 && (
          <Box sx={{ mb: 2 }} data-testid="delete-projects-with-files">
            <Typography variant="subtitle2">Folder with files: their files are kept.</Typography>
            {withFiles.map((project) => (
              <Typography key={project.name} variant="body2" sx={{ overflowWrap: 'anywhere' }}>
                {`${project.name} - ${project.filesDirectory}`}
              </Typography>
            ))}
          </Box>
        )}
        {noFolder.length > 0 && (
          <Box sx={{ mb: 2 }} data-testid="delete-projects-no-folder">
            <Typography variant="subtitle2">No folder or missing folder: nothing on disk to delete.</Typography>
            {noFolder.map((project) => (
              <Typography key={project.name} variant="body2">{project.name}</Typography>
            ))}
          </Box>
        )}
        <Typography variant="body2" sx={{ mb: 2 }}>
          projects.zip is saved first. If any export fails, nothing is deleted.
        </Typography>
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
        {step && <ExportProgress step={step} />}
        {deleting && (
          <Typography variant="body2" color="text.secondary" data-testid="delete-projects-deleting">
            Deleting the projects...
          </Typography>
        )}
      </DialogContent>
      <DialogActions>
        <Button onClick={onCancel} disabled={running}>
          Cancel
        </Button>
        <Button
          variant="contained"
          color="error"
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
