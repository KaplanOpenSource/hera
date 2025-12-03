import { Add } from "@mui/icons-material";
import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogContentText,
  DialogTitle,
  TextField,
} from "@mui/material";
import { useState } from "react";
import { execPython } from "../io/execPython";
import { fetchProjectDetails } from "../io/FetchProjects";
import { useProjectStore } from "../stores/useProjectStore";
import { ButtonTooltip } from "../elements/ButtonTooltip";
import { DocumentDesc, ProjectEntire, Toolkit } from "@shared/types";

export const AddDocumentButton = ({
  toolkit = undefined,
}: {
  toolkit?: Toolkit | undefined,
}) => {
  const [open, setOpen] = useState(false);
  const [name, setName] = useState('');
  const { currProjectName, setCurrentProject } = useProjectStore();

  const go = async () => {
    const desc: DocumentDesc = { datasourceName: name };
    if (toolkit?.toolkit) {
      desc.toolkit = toolkit.toolkit;
    }
    const { problem, data } = await execPython(`
import json
from hera.datalayer import All
All.addDocument('${currProjectName}', desc=${JSON.stringify(desc)})

docs = All.getDocumentsAsDict('${currProjectName}', with_id=True)
project = {"name": '${currProjectName}', "documents": docs['documents']}
result = json.dumps(project,indent=4)
      `)
    if (problem) {
      return;
    }

    const project = JSON.parse(data) as ProjectEntire;
    setCurrentProject(project);
  }

  return (<>
    <ButtonTooltip
      title='Add Document'
      onClick={() => setOpen(true)}
    >
      <Add />
    </ButtonTooltip>
    <Dialog open={open} onClose={() => setOpen(false)} onClick={e => e.stopPropagation()}>
      <DialogTitle>Add Document</DialogTitle>
      <DialogContent>
        <DialogContentText>
          Adding a document
        </DialogContentText>
        <TextField
          autoFocus
          required
          margin="dense"
          size="small"
          label="Name"
          fullWidth
          variant="outlined"
          value={name}
          onChange={(e) => setName(e.target.value)}
        />
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