import { Add } from "@mui/icons-material";
import {
  Button,
  DialogActions,
  DialogContent,
  DialogContentText,
  DialogTitle
} from "@mui/material";
import { DocumentDesc, ProjectEntire, Toolkit } from "@shared/types";
import { useRef, useState } from "react";
import { BooleanProperty } from "../../elements/BooleanProperty";
import { ButtonDialog } from "../../elements/ButtonDialog";
import { TextProperty } from "../../elements/TextProperty";
import { execPython } from "../../io/execPython";
import { useProjectStore } from "../../stores/useProjectStore";
import { SelectProperty } from "../../elements/SelectProperty";

const NO_TOOLKIT = '* No Toolkit *';

export const AddDocumentButton = ({
  toolkit = undefined,
}: {
  toolkit?: Toolkit | undefined,
}) => {
  const { toolkits } = useProjectStore();
  const [name, setName] = useState('');
  const [resource, setResource] = useState('');
  const [asAgent, setAsAgent] = useState(false);
  const [chosenToolkit, setChosenToolkit] = useState<string | undefined>(toolkit?.toolkit);
  const { currProjectName, setCurrentProject } = useProjectStore();
  const inputRef = useRef();
  const closeRef = useRef<() => void>();

  const doAddDoc = async () => {
    const desc: DocumentDesc = { datasourceName: name };
    if (chosenToolkit) {
      desc.toolkit = chosenToolkit;
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
            <TextProperty
              inputRef={inputRef}
              autoFocus
              required
              margin="dense"
              label="Name"
              fullWidth
              value={name}
              setValue={v => setName(v)}
            />
            <TextProperty
              margin="dense"
              label="Resource"
              fullWidth
              value={resource}
              setValue={v => setResource(v)}
              disabled={asAgent}
            />
            <BooleanProperty
              label="Agent"
              value={asAgent}
              setValue={v => setAsAgent(v)}
            />
            <SelectProperty
              label="Toolkit"
              value={chosenToolkit || NO_TOOLKIT}
              setValue={(v) => setChosenToolkit(v === NO_TOOLKIT ? undefined : v)}
              menuItems={[{ name: NO_TOOLKIT }, ...toolkits.map(t => ({ name: t.toolkit }))]}
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