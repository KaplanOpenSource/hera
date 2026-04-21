import { Add } from "@mui/icons-material";
import {
  Autocomplete,
  Button,
  DialogActions,
  DialogContent,
  DialogContentText,
  DialogTitle,
  Stack,
  TextField
} from "@mui/material";
import { DocumentDesc, ProjectEntire, Toolkit } from "../../shared/types";
import { useRef, useState } from "react";
import { ButtonDialog } from "../../elements/ButtonDialog";
import { TextProperty } from "../../elements/TextProperty";
import { fetchPython } from "../../io/fetchPython";
import { useProjectStore } from "../../stores/useProjectStore";
import { SelectProperty } from "../../elements/SelectProperty";

const METADATA_CLASSES = [
  { name: 'Metadata.Measurements', collection: 'Measurements_Collection' },
  { name: 'Metadata.Simulations', collection: 'Simulations_Collection' },
  { name: 'Metadata.Cache', collection: 'Cache_Collection' },
] as const;

type MetadataCls = typeof METADATA_CLASSES[number];

const NO_TOOLKIT = '* No Toolkit *';

const DOC_KINDS = ['Regular', 'Agent', 'Notebook'] as const;
type DocKind = typeof DOC_KINDS[number];

export const AddDocumentButton = ({
  toolkit = undefined,
  onDocumentCreated,
}: {
  toolkit?: Toolkit | undefined,
  onDocumentCreated?: (docOid: string) => void,
}) => {
  const { toolkits } = useProjectStore();
  const [name, setName] = useState('');
  const [resource, setResource] = useState('');
  const [kind, setKind] = useState<DocKind>('Regular');
  const [cls, setCls] = useState<MetadataCls>(METADATA_CLASSES[0]);
  const [chosenToolkit, setChosenToolkit] = useState<string | undefined>(toolkit?.toolkit);
  const { currProjectName, setCurrentProject } = useProjectStore();
  const inputRef = useRef();
  const closeRef = useRef<() => void>();

  const filesDir = useProjectStore.getState().getProject()?.configDocument?.data.desc.filesDirectory ?? '';

  const notebookResource = `${filesDir}/notebooks/${name}.ipynb`;

  const buildAddCommand = (desc: DocumentDesc) => {
    if (kind === 'Notebook') {
      return `
import json
from pathlib import Path
notebook_path = Path("${notebookResource}")
if not notebook_path.exists():
    notebook_path.parent.mkdir(parents=True, exist_ok=True)
    empty_notebook = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {"kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        }},
        "cells": [],
    }
    notebook_path.write_text(json.dumps(empty_notebook, indent=2))
Cache_Collection().addDocument(
    projectName="${currProjectName}",
    resource="${notebookResource}",
    dataFormat="JSON_dict",
    type="notebook",
    desc=${JSON.stringify(desc)},
)`;
    }
    if (kind === 'Agent') {
      return `
All.addDocument(
    '${currProjectName}',
    resource={"effects": {}},
    desc=${JSON.stringify(desc)},
    dataFormat=datatypes.JSON_DICT,
    type='ToolkitDataSource',
)`;
    }
    return `
${cls.collection}().addDocument(
    '${currProjectName}',
    resource='${resource}',
    desc=${JSON.stringify(desc)},
)`;
  };

  const doAddDoc = async () => {
    const desc: DocumentDesc = { datasourceName: name };
    if (chosenToolkit) {
      desc.toolkit = chosenToolkit;
    }
    const { data } = await fetchPython({
      results: ['project'],
      label: `add document ${name}`,
      code: `
from hera.datalayer import All, datatypes, Measurements_Collection, Simulations_Collection, Cache_Collection
${buildAddCommand(desc)}
docs = All.getDocumentsAsDict('${currProjectName}', with_id=True)
project = {"name": '${currProjectName}', "documents": docs['documents']}
`,
    })
    if (!data) {
      return;
    }
    const existingIds = useProjectStore.getState().getProject()?.documentIds ?? new Set<string>();
    setCurrentProject(data.project as ProjectEntire);
    const newIds = useProjectStore.getState().getProject()?.documentIds ?? new Set<string>();
    const [addedId] = newIds.difference(existingIds);
    if (addedId) {
      onDocumentCreated?.(addedId);
    }
  }

  return (
    <ButtonDialog
      icon={<Add />}
      title="Add Document"
      onOpen={() => {
        setName('');
        setResource('');
        setKind('Regular');
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
            <DialogContentText component="div" sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
              Adding a
              <SelectProperty
                label=""
                value={kind}
                setValue={v => {
                  setKind(v as DocKind);
                  if (v !== 'Regular') {
                    setChosenToolkit(undefined);
                  }
                }}
                menuItems={DOC_KINDS.map(k => ({ name: k }))}
              />
              document
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
              value={kind === 'Notebook' ? `${filesDir}/notebooks/${name}.ipynb` : resource}
              setValue={v => setResource(v)}
              disabled={kind !== 'Regular'}
              helperText={kind === 'Notebook' ? 'If a notebook file already exists at this path, it will be used as-is. Otherwise, a new empty notebook will be created.' : undefined}
            />
            {kind === 'Regular' && (
              <Stack
                direction={'column'}
                spacing={2}
                justifyItems={'flex-start'}
                alignItems={'flex-start'}
                mt={2}
              >
                <SelectProperty
                  label="Class"
                  value={cls.name}
                  setValue={(v) => setCls(METADATA_CLASSES.find(c => c.name === v)!)}
                  menuItems={METADATA_CLASSES.map(c => ({ name: c.name }))}
                />
                <Autocomplete
                  size="small"
                  disableClearable
                  style={{ minWidth: '200px' }}
                  options={[NO_TOOLKIT, ...toolkits.map(t => t.toolkit)]}
                  value={chosenToolkit || NO_TOOLKIT}
                  onChange={(_e, v) => setChosenToolkit(v === NO_TOOLKIT ? undefined : v)}
                  renderInput={(params) => <TextField {...params} label="Toolkit" />}
                  slotProps={{
                    popper: { placement: 'bottom-start', modifiers: [{ name: 'flip', enabled: false }] },
                    listbox: { style: { maxHeight: '250px' } },
                  }}
                />
              </Stack>
            )}
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