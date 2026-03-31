import { Add, Description } from '@mui/icons-material';
import { Stack, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { useCallback, useEffect, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { BASEURL } from '../../shared/baseurl';
import { idNotebookId, idNotebooksGroupId } from '../../shared/idDocId';

export const NotebookTreeItem = ({
  projectName,
  filesDir,
  onNotebookCreated,
}: {
  projectName: string,
  filesDir: string,
  onNotebookCreated: (itemId: string) => void,
}) => {
  const [notebooks, setNotebooks] = useState<string[]>([]);

  const refresh = useCallback(async () => {
    try {
      const r = await fetch(`${BASEURL}/notebooks?root_dir=${encodeURIComponent(filesDir)}`);
      const data = await r.json();
      setNotebooks(data.notebooks);
    } catch { /* ignore */ }
  }, [filesDir]);

  useEffect(() => { refresh(); }, [refresh]);

  const handleCreate = async () => {
    try {
      const r = await fetch(`${BASEURL}/notebooks/create`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ root_dir: filesDir }),
      });
      const data = await r.json();
      await refresh();
      onNotebookCreated(idNotebookId(data.name));
    } catch { /* ignore */ }
  };

  return (
    <TreeItem
      itemId={idNotebooksGroupId(projectName)}
      label={(
        <Stack direction="row" alignItems="center">
          <Typography marginRight={1}>Notebooks</Typography>
          <ButtonTooltip title="New Notebook" onClick={handleCreate}>
            <Add />
          </ButtonTooltip>
        </Stack>
      )}
      slots={{ icon: () => <Description fontSize="small" /> }}
    >
      {notebooks.map(name => (
        <TreeItem
          key={name}
          itemId={idNotebookId(name)}
          label={name}
          slots={{ icon: () => <Description fontSize="small" /> }}
        />
      ))}
    </TreeItem>
  );
};
