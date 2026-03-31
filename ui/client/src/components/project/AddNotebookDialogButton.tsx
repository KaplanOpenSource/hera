import { Add } from '@mui/icons-material';
import { Button, DialogActions, DialogContent, DialogTitle } from '@mui/material';
import { useRef, useState } from 'react';
import { ButtonDialog } from '../../elements/ButtonDialog';
import { TextProperty } from '../../elements/TextProperty';

export const AddNotebookDialogButton = ({
  defaultName,
  onCreate,
}: {
  defaultName: () => string,
  onCreate: (name: string) => void,
}) => {
  const [newName, setNewName] = useState('');
  const closeRef = useRef<() => void>();

  const doCreate = () => {
    const name = newName.trim();
    if (!name) return;
    closeRef.current?.();
    onCreate(name);
  };

  return (
    <ButtonDialog
      icon={<Add />}
      title="New Notebook"
      onOpen={() => setNewName(defaultName())}
      closeRef={closeRef}
      dialogProps={{
        onClick: (e) => e.stopPropagation(),
        onKeyDown: (e) => { if (e.code === 'Enter') doCreate(); },
      }}
    >
      {(close) => (<>
        <DialogTitle>New Notebook</DialogTitle>
        <DialogContent>
          <TextProperty
            autoFocus
            margin="dense"
            label="Name"
            fullWidth
            value={newName}
            setValue={setNewName}
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={close}>Cancel</Button>
          <Button onClick={doCreate} disabled={!newName.trim()}>Create</Button>
        </DialogActions>
      </>)}
    </ButtonDialog>
  );
};
