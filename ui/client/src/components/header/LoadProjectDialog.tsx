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
import { ProjectName } from "../../shared/types";
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
  const [readError, setReadError] = useState('');

  useEffect(() => {
    setZipInfo(null);
    setReadError('');
    file.arrayBuffer()
      .then((buffer) => setZipInfo(readProjectZip(new Uint8Array(buffer))))
      // A picked file may not be a hera export: say so in the dialog instead of failing silently.
      .catch((error) => {
        console.error(`could not read ${file.name}`, error);
        setReadError(`Could not read ${file.name} as a project export`);
      });
  }, [file]);

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
            value={zipInfo?.name ?? ''}
            fullWidth
            margin="dense"
            size="small"
            slotProps={{ input: { readOnly: true }, htmlInput: { 'data-testid': 'load-project-name' } }}
          />
          <TextField
            label="Files folder"
            value={zipInfo?.filesDirectory ?? ''}
            fullWidth
            margin="dense"
            size="small"
            slotProps={{ input: { readOnly: true }, htmlInput: { 'data-testid': 'load-project-folder' } }}
          />
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
          <Button onClick={() => loadProject()}>
            Load
          </Button>
        </DialogActions>
      </Dialog>
    </ThemeProvider>
  );
}
