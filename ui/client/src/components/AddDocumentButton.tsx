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
import { useState } from "react";
import { execPython } from "../io/execPython";
import { fetchProjectDetails, fetchProjectsNames } from "../io/FetchProjects";
import { useProjectStore } from "../stores/useProjectStore";
import { ButtonTooltip } from "../elements/ButtonTooltip";

export const AddDocumentButton = ({ }) => {
  const [open, setOpen] = useState(false);
  const [name, setName] = useState('');
  const [loadRepositories, setLoadRepositories] = useState(true);
  const { selectProject } = useProjectStore();

  const go = async () => {
    //     const { problem } = await execPython(`
    // from types import SimpleNamespace
    // from hera.utils.data.CLI import project_create
    // project_create(SimpleNamespace(projectName="${name}", directory=None, loadRepositories=${loadRepositories ? 'True' : 'False'}, overwrite=False))
    //       `)
    //     if (problem) {
    //       return;
    //     }
    //     await fetchProjectsNames();
    //     await fetchProjectDetails(name);
    //     selectProject(name);
  }

  return (<>
    <ButtonTooltip
      title='Add Document'
      onClick={() => setOpen(true)}
    >
      <Add />
    </ButtonTooltip>
    <Dialog open={open} onClose={() => setOpen(false)}>
      <DialogTitle>Add Document</DialogTitle>
      <DialogContent>
        <DialogContentText>
          Adding a document
        </DialogContentText>
        {/* <TextField
          autoFocus
          required
          margin="dense"
          size="small"
          label="Project Name"
          fullWidth
          variant="outlined"
          value={name}
          onChange={(e) => setName(e.target.value)}
        /> */}
        {/* <FormGroup>
          <FormControlLabel
            label="Load Repositories"
            control={<Checkbox
              checked={loadRepositories}
              onChange={(e) => setLoadRepositories(e.target.checked)}
            />}
          />
        </FormGroup> */}
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setOpen(false)}>
          Cancel
        </Button>
        <Button onClick={() => {
          go();
          setOpen(false);
        }}>
          Add Document
        </Button>
      </DialogActions>
    </Dialog>
  </>)
}