import { Folder, Settings } from "@mui/icons-material";
import { Stack, TextField, Tooltip, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";
import { repositoryApi } from "../../api/api";
import { CENTRAL_REPO_FOLDER_ID, idRepoId } from "../../shared/idDocId";

const DEFAULT_FOLDER = '~/hera/repositories/';
const STORAGE_KEY = 'hera-central-repo-folder';

const getFolder = () => localStorage.getItem(STORAGE_KEY) || DEFAULT_FOLDER;

export const CentralRepoFolder = () => {
  const [folder, setFolder] = useState(getFolder);
  const [files, setFiles] = useState<string[]>([]);
  const [loaded, setLoaded] = useState(false);
  const { openDialog, DialogComponent } = useDialog<{ folder: string }>();

  const fetchFiles = async (folderPath: string) => {
    const answer = await repositoryApi.jsonsInFolder(folderPath);
    if (answer) {
      setFiles(answer.files);
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
        itemId={CENTRAL_REPO_FOLDER_ID}
        label={
          <Tooltip title="Central repository folder. Contains local JSON repository files on disk that can be loaded into any project. Click the chevron to expand.">
            <Stack
              direction="row"
              spacing={1}
              alignItems="center"
            >
              <span>
                <b>Central folder</b>
              </span>
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
          ? <TreeItem itemId={`${CENTRAL_REPO_FOLDER_ID}/__loading`} label="Loading..." />
          : files.length === 0
            ? <TreeItem itemId={`${CENTRAL_REPO_FOLDER_ID}/__empty`} label="No JSON files found" />
            : files.map(f => (
              <TreeItem key={idRepoId(f)} itemId={idRepoId(f)} label={f} />
            ))
        }
      </TreeItem>
      {DialogComponent}
    </>
  )
}
