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
import { BooleanProperty } from "../../elements/BooleanProperty";
import { ButtonDialog } from "../../elements/ButtonDialog";

export const AddDocumentButton = ({
  toolkit = undefined,
}: {
  toolkit?: Toolkit | undefined,
}) => {
  const [open, setOpen] = useState(false);
  const [name, setName] = useState('');
  const [resource, setResource] = useState('');
  const [asAgent, setAsAgent] = useState(false);
  const { currProjectName, setCurrentProject } = useProjectStore();
  const inputRef = useRef();
  const closeRef = useRef<() => void>();

  const doAddDoc = async () => {
    const desc: DocumentDesc = { datasourceName: name };
    if (toolkit?.toolkit) {
      desc.toolkit = toolkit.toolkit;
    }
    const addcmd = asAgent
      ? `
All.addDocument('${currProjectName}', resource={"effects": {}}, desc=${JSON.stringify(desc)},
    dataFormat=datatypes.JSON_DICT,
    type='ToolkitDataSource')
      `
      : `
All.addDocument('${currProjectName}', resource='${resource}', desc=${JSON.stringify(desc)})
    `;
    const { problem, data } = await execPython(`
import json
from hera.datalayer import All, datatypes
${addcmd}

docs = All.getDocumentsAsDict('${currProjectName}', with_id=True)
result = {"name": '${currProjectName}', "documents": docs['documents']}
      `)
    if (problem) {
      return;
    }
    setCurrentProject((data as ProjectEntire));
  }

  return (
    <ButtonDialog
      icon={<Add />}
      title="Add Document"
      onOpen={() => {
        setName('');
        setResource('');
        setAsAgent(false);
        setTimeout(() => (inputRef.current as any)?.focus(), 0)

      }}
      dialogProps={{
        onClick: (e) => e.stopPropagation(),
        onKeyDown: (e) => {
          if (e.code === 'Enter') {
            doAddDoc();
            closeRef.current?.();
          }
        },
      }}
      closeRef={closeRef}
    >
      {(close) => (
        <>
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
              disabled={asAgent}
            />
            <BooleanProperty
              label="Agent"
              value={asAgent}
              setValue={v => setAsAgent(v)}
            />
          </DialogContent>
          <DialogActions>
            <Button onClick={close}>
              Cancel
            </Button>
            <Button onClick={() => {
              doAddDoc();
              close();
            }}>
              Add Document
            </Button>
          </DialogActions>
        </>
      )}
    </ButtonDialog>
  )
}