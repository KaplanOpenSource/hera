import { Add } from "@mui/icons-material";
import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogContentText,
  DialogTitle,
  TextField,
  ThemeProvider
} from "@mui/material";
import { useRef, useState } from "react";
import { useNavigate } from "react-router-dom";
import { BooleanProperty } from "../../elements/BooleanProperty";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { postJson } from "../../io/entrypoint";
import { SimpleTreeView } from "@mui/x-tree-view";
import { RegisteredRepositories } from "../repo/RegisteredRepositories";
import { useAppTheme } from "../../theme";
import { useProjectListStore } from "../../stores/useProjectListStore";
import { ProjectName } from "../../shared/types";

export const AddProjectButton = ({ }) => {
  // Follow the app theme so the dialog isn't tinted by the dark header it opens from.
  const dialogTheme = useAppTheme();
  const [open, setOpen] = useState(false);
  const [name, setName] = useState('');
  const [filesDirectory, setFilesDirectory] = useState('');
  const [loadRepositories, setLoadRepositories] = useState(true);
  const navigate = useNavigate();
  const inputRef = useRef();

  const trimmedName = name.trim();

  const doAddProject = async () => {
    if (!trimmedName) return;

    // The server creates the project and registers the notebooks already there.
    const data = await postJson('/project/create', {
      projectName: trimmedName,
      filesDirectory,
      loadRepositories,
    }, { label: `create project ${trimmedName}` });
    if (!data) {
      return;
    }

    // The URL stays the only thing that says which project is current, so the
    // project travels with the navigation: the watcher shows it from there
    // instead of loading it again. See FetchProjects.
    useProjectListStore.getState().setProjectNames(data.projects as ProjectName[]);
    navigate('/' + encodeURIComponent(trimmedName), { state: { project: data.project } });
    setOpen(false);
  }

  return (<>
    <ButtonTooltip
      title='Add project'
      onClick={() => {
        setOpen(true);
        setName('');
        setTimeout(() => (inputRef.current as any)?.focus(), 0)
      }}
    >
      <Add />
    </ButtonTooltip>
    <ThemeProvider theme={dialogTheme}>
    <Dialog
      open={open}
      onClose={() => setOpen(false)}
      onKeyDown={e => { if (e.code === 'Enter' && trimmedName) doAddProject() }}
    >
      <DialogTitle>New Project</DialogTitle>
      <DialogContent>
        <DialogContentText>
          Adding a new project, please fill initial details
        </DialogContentText>
        <TextField
          inputRef={inputRef}
          // autoFocus
          required
          margin="dense"
          size="small"
          label="Project Name"
          fullWidth
          variant="outlined"
          value={name}
          onChange={(e) => setName(e.target.value)}
        />
        <TextField
          margin="dense"
          size="small"
          label={"Project Files Directory (leave empty for default)"}
          fullWidth
          variant="outlined"
          value={filesDirectory}
          onChange={(e) => setFilesDirectory(e.target.value)}
        />
        <BooleanProperty
          label="Load Repositories"
          value={loadRepositories}
          setValue={v => setLoadRepositories(v)}
        />
        <SimpleTreeView defaultExpandedItems={['registered-repos']}>
          <RegisteredRepositories />
        </SimpleTreeView>
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setOpen(false)}>
          Cancel
        </Button>
        <Button onClick={() => doAddProject()} disabled={!trimmedName}>
          Add Project
        </Button>
      </DialogActions>
    </Dialog>
    </ThemeProvider>
  </>)
}