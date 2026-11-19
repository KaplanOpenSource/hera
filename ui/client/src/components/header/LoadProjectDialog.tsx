import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Stack,
  TextField,
  ThemeProvider,
  Typography
} from "@mui/material";
import { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import { projectApi } from "../../api/api";
import { useProjectListStore } from "../../stores/useProjectListStore";
import { ProjectFolder, ProjectName } from "../../shared/types";
import { useAppTheme } from "../../theme";
import { ProjectZipInfo, readProjectZip } from "../../utils/readProjectZip";

// Shows the picked zip, and loads it as a project.
export const LoadProjectDialog = ({
  file,
  onPickAnother,
  onClose,
  onLoaded,
}: {
  file: File,
  onPickAnother: () => void,
  onClose: () => void,
  onLoaded?: () => void,
}) => {
  // Follow the app theme so the dialog isn't tinted by the dark header it opens from.
  const dialogTheme = useAppTheme();
  const navigate = useNavigate();
  const [zipInfo, setZipInfo] = useState<ProjectZipInfo | null>(null);
  // Starts as what the zip holds; edits only change the checks, not what Load sends yet.
  const [edited, setEdited] = useState<ProjectZipInfo>({ name: '', filesDirectory: '' });
  const [readError, setReadError] = useState('');
  const [folders, setFolders] = useState<ProjectFolder[]>([]);

  useEffect(() => {
    projectApi.folders().then((data) => setFolders(data?.projects ?? []));
  }, []);

  useEffect(() => {
    setZipInfo(null);
    setEdited({ name: '', filesDirectory: '' });
    setReadError('');
    file.arrayBuffer()
      .then((buffer) => {
        const info = readProjectZip(new Uint8Array(buffer));
        setZipInfo(info);
        setEdited(info);
      })
      // A picked file may not be a hera export: say so in the dialog instead of failing silently.
      .catch((error) => {
        console.error(`could not read ${file.name}`, error);
        setReadError(`Could not read ${file.name} as a project export`);
      });
  }, [file]);

  const nameTaken = !!zipInfo && folders.some((p) => p.name === edited.name);
  const sharedWith = folders
    .filter((p) => !!zipInfo && !!p.filesDirectory && p.filesDirectory === edited.filesDirectory)
    .map((p) => p.name);

  const loadProject = async () => {
    const data = await projectApi.load(await file.arrayBuffer(), file.name);
    if (!data) {
      return;
    }
    // Same as after create: the project travels with the navigation. See AddProjectButton.
    useProjectListStore.getState().setProjectNames(data.projects as ProjectName[]);
    navigate('/' + encodeURIComponent(data.project.name), { state: { project: data.project } });
    onClose();
    onLoaded?.();
  }

  return (
    <ThemeProvider theme={dialogTheme}>
      <Dialog
        open
        onClose={onClose}
      >
        <DialogTitle>Load a project</DialogTitle>
        <DialogContent>
          <Stack direction="row" alignItems="center" spacing={1}>
            <Typography data-testid="load-project-file-name">
              {file.name}
            </Typography>
            <Button
              variant="outlined"
              onClick={onPickAnother}
            >
              Choose another zip
            </Button>
          </Stack>
          <TextField
            label="Project name"
            value={edited.name}
            onChange={(e) => setEdited({ ...edited, name: e.target.value })}
            fullWidth
            margin="dense"
            size="small"
            slotProps={{ input: { readOnly: !zipInfo }, htmlInput: { 'data-testid': 'load-project-name' } }}
          />
          {nameTaken && (
            <Typography variant="caption" color="error" data-testid="load-project-name-taken">
              A project named {edited.name} already exists
            </Typography>
          )}
          <TextField
            label="Files folder"
            value={edited.filesDirectory}
            onChange={(e) => setEdited({ ...edited, filesDirectory: e.target.value })}
            fullWidth
            margin="dense"
            size="small"
            slotProps={{ input: { readOnly: !zipInfo }, htmlInput: { 'data-testid': 'load-project-folder' } }}
          />
          {sharedWith.length > 0 && (
            <Typography variant="caption" color="warning.main" data-testid="load-project-shared-folder">
              This folder is also used by: {sharedWith.join(', ')}
            </Typography>
          )}
          {readError && (
            <Typography color="error" data-testid="load-project-read-error">
              {readError}
            </Typography>
          )}
        </DialogContent>
        <DialogActions>
          <Button onClick={onClose}>
            Cancel
          </Button>
          <Button onClick={() => loadProject()} disabled={!zipInfo || nameTaken}>
            Load
          </Button>
        </DialogActions>
      </Dialog>
    </ThemeProvider>
  );
}
