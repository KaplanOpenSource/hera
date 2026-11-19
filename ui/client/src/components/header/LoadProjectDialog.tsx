import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Stack,
  ThemeProvider,
  Typography
} from "@mui/material";
import { useNavigate } from "react-router-dom";
import { projectApi } from "../../api/api";
import { useProjectListStore } from "../../stores/useProjectListStore";
import { ProjectName } from "../../shared/types";
import { useAppTheme } from "../../theme";

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
