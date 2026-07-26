import { Box, Stack, Typography } from '@mui/material';
import { useState } from 'react';
import { RenameField } from '../../elements/RenameField';

export const DetailsViewItemName = ({
  itemKey,
  setItemKey = undefined,
}: {
  itemKey: string,
  setItemKey?: (newKey: string | undefined) => void | undefined,
}) => {
  const [editing, setEditing] = useState(false);

  return (
    <Stack direction={'row'} alignItems={'center'}>
      {(editing && setItemKey)
        ? (
          <Box
            sx={{
              minWidth: '300px'
            }}
          >
            <RenameField
              value={itemKey}
              setValue={val => {
                setItemKey(val);
                setEditing(false);
              }}
            />
          </Box>
        )
        : (
          <Typography
            onClick={() => setItemKey && setEditing(true)}
            sx={{
              whiteSpace: 'nowrap',
              minWidth: '100px',
              cursor: setItemKey ? 'text' : 'default'
            }}
          >
            {itemKey}
          </Typography>
        )}
    </Stack>
  )
}
