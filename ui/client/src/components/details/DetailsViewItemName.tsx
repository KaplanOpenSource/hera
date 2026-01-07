import { Edit } from '@mui/icons-material';
import { Box, Stack, Typography } from '@mui/material';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
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
        : (<>
          <Typography
            sx={{
              whiteSpace: 'nowrap',
              minWidth: '100px'
            }}
          >
            {itemKey}
          </Typography>
          {setItemKey && (<>
            <ButtonTooltip
              title={'Rename ' + itemKey}
              onClick={() => setEditing(!editing)}
            >
              <Edit />
            </ButtonTooltip>
          </>)}
        </>)}
    </Stack>
  )
}
