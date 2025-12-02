import { Add } from "@mui/icons-material";
import {
  Button,
  Checkbox,
  Dialog,
  DialogActions,
  DialogContent,
  DialogContentText,
  DialogTitle,
  FormControlLabel,
  FormGroup,
  TextField,
} from "@mui/material";
import { useRef, useState } from "react";
import { execPython } from "../io/execPython";
import { fetchProjectDetails, fetchProjectsNames } from "../io/FetchProjects";
import { useProjectStore } from "../stores/useProjectStore";
import { ButtonTooltip } from "../elements/ButtonTooltip";

export const AddProjectButton = ({ }) => {
  const [open, setOpen] = useState(false);
  const [name, setName] = useState('');
  const [loadRepositories, setLoadRepositories] = useState(true);
  const { selectProject, projectNames } = useProjectStore();
  const inputRef = useRef();

  const doAddProject = async () => {
    const { problem } = await execPython(`
import os
from types import SimpleNamespace
from hera.utils.data.CLI import project_create
project_create(SimpleNamespace(
  projectName='${name}',
  directory=os.path.join(os.getcwd(), 'projects', '${name}'),
  loadRepositories=${loadRepositories ? 'True' : 'False'},
  overwrite=False))
      `)
    if (problem) {
      return;
    }
    setTimeout(async () => {
      // debugger
      await fetchProjectsNames();
      // debugger
      setTimeout(async () => {
        console.log('after fetch', projectNames.map(x => x.name))
        await fetchProjectDetails(name);
        setTimeout(async () => {
          selectProject(name);
        }, 0);
      }, 0);
    }, 0);
    setOpen(false);
  }

  return (<>
    <ButtonTooltip
      title='Add project'
      onClick={() => {
        setOpen(true);
        setTimeout(() => (inputRef.current as any)?.focus(), 0)
      }}
    >
      <Add />
    </ButtonTooltip>
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
        <FormGroup>
          <FormControlLabel
            label="Load Repositories"
            control={<Checkbox
              checked={loadRepositories}
              onChange={(e) => setLoadRepositories(e.target.checked)}
            />}
          />
        </FormGroup>
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
  </>)
}