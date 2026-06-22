import { Add } from "@mui/icons-material";
import {
  Autocomplete,
  Button,
  DialogActions,
  DialogContent,
  DialogTitle,
  Stack,
  TextField,
  ToggleButton,
  ToggleButtonGroup,
} from "@mui/material";
import { DocumentDesc, ProjectEntire, Toolkit } from "../../shared/types";
import { useRef, useState } from "react";
import { ButtonDialog } from "../../elements/ButtonDialog";
import { TextProperty } from "../../elements/TextProperty";
import { fetchPython } from "../../io/fetchPython";
import { useProjectStore } from "../../stores/useProjectStore";
import { SelectProperty } from "../../elements/SelectProperty";
import { buildAddDocumentCode } from "./buildAddDocumentCode";

const METADATA_CLASSES = [
  { name: 'Metadata.Measurements', collection: 'Measurements_Collection' },
  { name: 'Metadata.Simulations', collection: 'Simulations_Collection' },
  { name: 'Metadata.Cache', collection: 'Cache_Collection' },
] as const;

type MetadataCls = typeof METADATA_CLASSES[number];

const NO_TOOLKIT = '* No Toolkit *';

export enum DocKind {
  Document = 'Document',
  Agent = 'Agent',
  Notebook = 'Notebook',
  Workflow = 'Workflow',
}

const DOC_KINDS = Object.values(DocKind);

const getNextDefaultName = (kind: DocKind) => {
  const prefix = kind;
  const regex = new RegExp(`^${prefix}(\\d+)$`);
  const docs = useProjectStore.getState().getProject()?.documents ?? [];
  let max = 0;
  for (const doc of docs) {
    const match = doc.name.match(regex);
    if (match) {
      max = Math.max(max, parseInt(match[1], 10));
    }
  }
  return `${prefix}${max + 1}`;
};

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
  const [kind, setKind] = useState(DocKind.Document);
  const [cls, setCls] = useState<MetadataCls>(METADATA_CLASSES[0]);
  const [chosenToolkit, setChosenToolkit] = useState<string | undefined>(toolkit?.toolkit);
  const { currProjectName, setCurrentProject } = useProjectStore();
  const inputRef = useRef();
  const closeRef = useRef<() => void>();

  const filesDir = useProjectStore.getState().getProject()?.configDocument?.data.desc.filesDirectory ?? '';

  const notebookResource = `${filesDir}/notebooks/${name}.ipynb`;

  const doAddDoc = async () => {
    const desc: DocumentDesc = { datasourceName: name };
    if (chosenToolkit) {
      desc.toolkit = chosenToolkit;
    }
    const { data } = await fetchPython({
      results: ['project'],
      label: `add document ${name}`,
      code: buildAddDocumentCode({
        kind,
        projectName: currProjectName,
        desc,
        toolkitNames: useProjectStore.getState().getProjectToolkitKeys(),
        notebookResource,
        collection: cls.collection,
        resource,
      }),
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
      icon={
        <Stack direction={'row'} alignItems={'center'} spacing={0.5}>
          <Add fontSize="small" />
        </Stack>
      }
      button
      title="Add document"
      onOpen={() => {
        setName(getNextDefaultName(DocKind.Document));
        setResource('');
        setKind(DocKind.Document);
        setTimeout(() => (inputRef.current as any)?.select(), 0)
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
            <ToggleButtonGroup
              exclusive
              size="small"
              value={kind}
              onChange={(_e, v) => {
                if (v === null) return;
                const newKind = v as DocKind;
                setKind(newKind);
                setName(getNextDefaultName(newKind));
                if (newKind !== DocKind.Document) {
                  setChosenToolkit(undefined);
                }
                setTimeout(() => (inputRef.current as any)?.select(), 0);
              }}
              sx={{ mb: 1 }}
            >
              {DOC_KINDS.map(k => (
                <ToggleButton key={k} value={k}>{k}</ToggleButton>
              ))}
            </ToggleButtonGroup>
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
              value={kind === DocKind.Notebook ? `${filesDir}/notebooks/${name}.ipynb` : resource}
              setValue={v => setResource(v)}
              disabled={kind !== DocKind.Document}
              helperText={kind === DocKind.Notebook ? 'If a notebook file already exists at this path, it will be used as-is. Otherwise, a new empty notebook will be created.' : undefined}
            />
            {kind === DocKind.Document && (
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