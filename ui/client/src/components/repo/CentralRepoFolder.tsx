import { Folder, Settings } from "@mui/icons-material";
import { Stack, TextField, Tooltip, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";
import { fetchPython } from "../../io/fetchPython";
import { idRepoId } from "../../shared/idDocId";
import { LOAD_REPO_JSON_PYTHON } from "../../shared/repoJsonPython";

const DEFAULT_FOLDER = '~/hera/repositories/';
const STORAGE_KEY = 'hera-central-repo-folder';

const getFolder = () => localStorage.getItem(STORAGE_KEY) || DEFAULT_FOLDER;

export const CentralRepoFolder = () => {
  const [folder, setFolder] = useState(getFolder);
  const [files, setFiles] = useState<string[]>([]);
  const [loaded, setLoaded] = useState(false);
  const { openDialog, DialogComponent } = useDialog<{ folder: string }>();

  const fetchFiles = async (folderPath: string) => {
    const { data } = await fetchPython({
      results: ['jsonFiles'],
      label: 'central repo files',
      code: `
import os, glob, json
${LOAD_REPO_JSON_PYTHON}
folder = os.path.expanduser('${folderPath}')
jsonFiles = []
if os.path.isdir(folder):
    allFiles = glob.glob(os.path.join(folder, '**', '*.json'), recursive=True)
    allFiles = [f for f in allFiles if not f.endswith('caseConfiguration.json')]
    jsonFiles = sorted(f for f in allFiles if loadRepoJson(f) is not None)
`,
    });
    if (data?.jsonFiles) {
      setFiles(data.jsonFiles);
    }
    setLoaded(true);
  };

  const handleExpand = async () => {
    if (loaded) return;
    await fetchFiles(folder);
  };

  const handleChangeFolder = async () => {
    const result = await openDialog({
      title: 'Central Repository Folder',
      initialValues: { folder },
      yesText: 'Save',
      noText: 'Cancel',
      render: ({ values, setValues }) => (
        <TextField
          autoFocus
          fullWidth
          label="Folder path"
          value={values.folder}
          onChange={(e) => setValues({ folder: e.target.value })}
          onClick={(e) => e.stopPropagation()}
          onKeyDown={(e) => e.stopPropagation()}
          size="small"
        />
      ),
    });
    if (result.confirmed && result.values?.folder) {
      localStorage.setItem(STORAGE_KEY, result.values.folder);
      setFolder(result.values.folder);
      setLoaded(false);
      setFiles([]);
      await fetchFiles(result.values.folder);
    }
  };

  return (
    <>
      <TreeItem
        itemId="central-repo-folder"
        label={
          <Tooltip title="Central repository folder. Contains local JSON repository files on disk that can be loaded into any project. Click the chevron to expand.">
            <Stack
              direction="row"
              spacing={1}
              alignItems="center"
            >
              <Folder fontSize="small" />
              <Typography variant="body2">{folder}</Typography>
              <ButtonTooltip title="Change repository folder"
                onClick={(e) => {
                  e.stopPropagation();
                  handleChangeFolder();
                }}
              >
                <Settings fontSize="small" />
              </ButtonTooltip>
            </Stack>
          </Tooltip>
        }
        onClick={handleExpand}
      >
        {!loaded
          ? <TreeItem itemId="central-repo-folder/__loading" label="Loading..." />
          : files.length === 0
            ? <TreeItem itemId="central-repo-folder/__empty" label="No JSON files found" />
            : files.map(f => (
              <TreeItem key={idRepoId(f)} itemId={idRepoId(f)} label={f} />
            ))
        }
      </TreeItem>
      {DialogComponent}
    </>
  )
}
