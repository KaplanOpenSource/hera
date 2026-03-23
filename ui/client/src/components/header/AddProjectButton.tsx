import { Add } from "@mui/icons-material";
import {
  Button,
  createTheme,
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
import { fetchPython } from "../../io/fetchPython";
import { RepositoriesInProject } from "../repo/RepositoriesInProject";

export const AddProjectButton = ({ }) => {
  const [open, setOpen] = useState(false);
  const [name, setName] = useState('');
  const [filesDirectory, setFilesDirectory] = useState('');
  const [loadRepositories, setLoadRepositories] = useState(true);
  const navigate = useNavigate();
  const inputRef = useRef();

  const doAddProject = async () => {
    let dirStr = `os.path.join(os.getcwd(), 'projects', '${name}')`;
    if (filesDirectory !== '') {
      dirStr = `'${filesDirectory}'`;
    }

    const { problem } = await fetchPython({
      results: [],
      code: `
import os
from types import SimpleNamespace

from hera.utils.data.CLI import project_create
from hera.datalayer.project import Project

project_create(SimpleNamespace(
  projectName='${name}',
  directory=${dirStr},
  loadRepositories=${loadRepositories ? 'True' : 'False'},
  overwrite=False))

# Creates a config document in MongoDB so the project appears in getProjectList()
Project(projectName='${name}', filesDirectory=${dirStr})
`,
    })
    if (problem) {
      return;
    }

    navigate('/' + encodeURIComponent(name));
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
    <ThemeProvider theme={createTheme()}>
    <Dialog
      open={open}
      onClose={() => setOpen(false)}
      onKeyDown={e => { if (e.code === 'Enter') doAddProject() }}
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
        <RepositoriesInProject />
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setOpen(false)}>
          Cancel
        </Button>
        <Button onClick={() => doAddProject()}>
          Add Project
        </Button>
      </DialogActions>
    </Dialog>
    </ThemeProvider>
  </>)
}