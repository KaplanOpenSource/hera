import { Add } from "@mui/icons-material";
import { Button, Dialog, DialogActions, DialogContent, DialogContentText, DialogTitle, IconButton, TextField } from "@mui/material";
import { useState } from "react";
import { execPython } from "../io/execPython";
import { fetchProjectDetails, fetchProjectsNames } from "../io/FetchProjects";
import { useProjectStore } from "../stores/useProjectStore";

export const AddProjectButton = ({ }) => {
  const [open, setOpen] = useState(false);
  const [name, setName] = useState('');
  const { selectProject } = useProjectStore();

  const addProject = async () => {
    const { problem } = await execPython(`
from types import SimpleNamespace
from hera.utils.data.CLI import project_create
project_create(SimpleNamespace(projectName="${name}", directory=None, loadRepositories=True, overwrite=False))
      `)
    if (problem) {
      return;
    }
    await fetchProjectsNames();
    await fetchProjectDetails(name);
    selectProject(name);
  }

  return (<>
    <IconButton
      onClick={() => setOpen(true)}
    >
      <Add />
    </IconButton>
    <Dialog open={open} onClose={() => setOpen(false)}>
      <DialogTitle>New Project</DialogTitle>
      <DialogContent>
        <DialogContentText>
          Adding a new project, please fill initial details
        </DialogContentText>
        <TextField
          autoFocus
          required
          margin="dense"
          size="small"
          label="Project Name"
          fullWidth
          variant="outlined"
          value={name}
          onChange={(e) => setName(e.target.value)}
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setOpen(false)}>
          Cancel
        </Button>
        <Button onClick={() => {
          addProject();
          setOpen(false);
        }}>
          Add Project
        </Button>
      </DialogActions>
    </Dialog>
  </>)
}