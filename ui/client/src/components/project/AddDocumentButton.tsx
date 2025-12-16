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
import { useRef, useState } from "react";
import { execPython } from "../../io/execPython";
import { useProjectStore } from "../../stores/useProjectStore";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { DocumentDesc, ProjectEntire, Toolkit } from "@shared/types";

export const AddDocumentButton = ({
  toolkit = undefined,
}: {
  toolkit?: Toolkit | undefined,
}) => {
  const [open, setOpen] = useState(false);
  const [name, setName] = useState('');
  const [resource, setResource] = useState('');
  const { currProjectName, setCurrentProject } = useProjectStore();
  const inputRef = useRef();

  const doAddDoc = async () => {
    const desc: DocumentDesc = { datasourceName: name };
    if (toolkit?.toolkit) {
      desc.toolkit = toolkit.toolkit;
    }
    const { problem, data } = await execPython(`
import json
from hera.datalayer import All
All.addDocument('${currProjectName}', resource='${resource}', desc=${JSON.stringify(desc)})

docs = All.getDocumentsAsDict('${currProjectName}', with_id=True)
result = {"name": '${currProjectName}', "documents": docs['documents']}
      `)
    if (problem) {
      return;
    }
    setCurrentProject((data as ProjectEntire));
    setOpen(false);
  }

  return (<>
    <ButtonTooltip
      title='Add Document'
      onClick={() => {
        setOpen(true)
        setName('');
        setTimeout(() => (inputRef.current as any)?.focus(), 0)
      }}
    >
      <Add />
    </ButtonTooltip>
    <Dialog
      open={open}
      onClose={() => setOpen(false)}
      onClick={e => e.stopPropagation()}
      onKeyDown={e => { if (e.code === 'Enter') doAddDoc() }}
    >
      <DialogTitle>Add Document</DialogTitle>
      <DialogContent>
        <DialogContentText>
          Adding a document
        </DialogContentText>
        <TextField
          inputRef={inputRef}
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
        <TextField
          margin="dense"
          size="small"
          label="Resource"
          fullWidth
          variant="outlined"
          value={resource}
          onChange={(e) => setResource(e.target.value)}
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setOpen(false)}>
          Cancel
        </Button>
        <Button onClick={() => {
          doAddDoc();
        }}>
          Add Document
        </Button>
      </DialogActions>
    </Dialog>
  </>)
}