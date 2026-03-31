import { Add, Refresh } from '@mui/icons-material';
import { Stack, Typography } from '@mui/material';
import { TreeItem } from '@mui/x-tree-view';
import { useCallback, useEffect, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { fetchPython } from '../../io/fetchPython';
import { NotebookCommands } from '../../io/NotebookCommands';
import { ID_NOTEBOOKS_GROUP, idFromNotebookId, idNotebookId } from '../../shared/idDocId';
import { NotebookListItem } from './NotebookListItem';

export const NotebookList = ({
  filesDir,
  selectedItemId,
  setSelectedItemIds,
}: {
  filesDir: string,
  selectedItemId: string | undefined,
  setSelectedItemIds: (ids: string[]) => void,
}) => {
  const [notebooks, setNotebooks] = useState<string[]>([]);

  const refresh = useCallback(async () => {
    const { data } = await fetchPython(NotebookCommands.list(filesDir));
    if (data) setNotebooks(data.notebooks);
  }, [filesDir]);

  useEffect(() => { refresh(); }, [refresh]);

  const handleCreate = async () => {
    const { data } = await fetchPython(NotebookCommands.create(filesDir));
    await refresh();
    if (data) setSelectedItemIds([idNotebookId(data.name)]);
  };

  const handleDelete = async (name: string) => {
    await fetchPython(NotebookCommands.delete(filesDir, name));
    await refresh();
    if (selectedItemId === idNotebookId(name)) {
      setSelectedItemIds([]);
    }
  };

  const handleItemRefresh = useCallback(() => {
    if (!selectedItemId) return;
    const name = idFromNotebookId(selectedItemId);
    if (name && notebooks.includes(name)) {
      setSelectedItemIds([]);
      setTimeout(() => setSelectedItemIds([selectedItemId]), 0);
    } else {
      setSelectedItemIds([]);
    }
  }, [selectedItemId, notebooks, setSelectedItemIds]);

  return (
    <TreeItem
      itemId={ID_NOTEBOOKS_GROUP}
      label={(
        <Stack direction="row" alignItems="center">
          <Typography marginRight={1}>Notebooks</Typography>
          <ButtonTooltip title="Refresh" onClick={refresh}>
            <Refresh />
          </ButtonTooltip>
          <ButtonTooltip title="New Notebook" onClick={handleCreate}>
            <Add />
          </ButtonTooltip>
        </Stack>
      )}
    >
      {notebooks.map(name => (
        <NotebookListItem
          key={name}
          name={name}
          selectedItemId={selectedItemId}
          onDelete={handleDelete}
          onRefresh={handleItemRefresh}
        />
      ))}
    </TreeItem>
  );
};
