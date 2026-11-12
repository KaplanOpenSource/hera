import { Box, Typography } from '@mui/material';
import { ReactNode } from 'react';
import { RenameField } from '../../elements/RenameField';
import { DeleteFieldButton } from './DeleteFieldButton';
import { DESC_FIELD } from '../../shared/constants';

// A row's name, with its delete button, in the details tree. The name is the
// raw key unless `nameForView` gives a friendlier one to show.
export const DetailsViewItemName = ({
  itemKey,
  parentKey,
  setItemKey = undefined,
  nameForView = undefined,
  allowRename = true,
}: {
  itemKey: string,
  parentKey?: string,
  setItemKey?: (newKey: string | undefined) => void | undefined,
  nameForView?: (itemKey: string, parentKey: string | undefined) => ReactNode,
  allowRename?: boolean,
}) => {
  // The top-level `desc` field isn't renameable, so it gets a label of its own.
  let shownName = nameForView?.(itemKey, parentKey);
  if (shownName === undefined && itemKey === DESC_FIELD && !parentKey) {
    shownName = (
      <Typography sx={{ whiteSpace: 'nowrap', minWidth: '100px', flexShrink: 0 }}>
        Description (desc)
      </Typography>
    );
  }

  return (
    // The delete button sits on the name's top-left corner, over no text.
    <Box
      sx={{
        position: 'relative',
        display: 'flex',
        minWidth: 0,
        // A row with a custom name keeps it whole; other names may be cut short.
        flexShrink: shownName !== undefined ? 0 : 1,
      }}
    >
      {!allowRename && shownName}
      {allowRename && (
        <RenameField
          value={itemKey}
          setValue={setItemKey}
          labelMinWidth="100px"
          valueForView={shownName}
          // A custom name is shown in full; the value field gives up the room.
          keepViewWidth
        />
      )}
      {setItemKey && (
        <Box
          className="field-delete"
          sx={{
            position: 'absolute',
            left: '-8px',
            top: 0,
            transform: 'translateY(-40%)',
            zIndex: 2,
            '& .MuiIconButton-root': { padding: '2px' },
            '& .MuiSvgIcon-root': { fontSize: '0.8rem' },
          }}
        >
          <DeleteFieldButton
            itemKey={itemKey}
            setItemKey={setItemKey}
          />
        </Box>
      )}
    </Box>
  );
};
