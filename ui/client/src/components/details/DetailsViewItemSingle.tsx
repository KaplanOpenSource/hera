import { FormControl, InputLabel, MenuItem, Select, SelectChangeEvent, Stack, TextField, Typography } from '@mui/material';
import { Case, SwitchCase } from '../../elements/SwitchCase';
import { useServerConstants } from '../../stores/useServerConstants';

export const DetailsViewItemSingle = ({
  itemKey,
  itemValue,
  setItemValue,
}: {
  itemKey: string,
  itemValue: any,
  setItemValue: (newVal: string) => void,
}) => {
  const { dataTypes } = useServerConstants();
  return (
    <Stack direction='row' spacing={1} justifyItems={'center'} alignItems={'center'}>
      <Typography>
        {itemKey}
      </Typography>
      <SwitchCase test={itemKey}>
        <Case isDefault>
          <TextField
            size='small'
            value={['object', 'function'].includes(typeof itemValue) ? JSON.stringify(itemValue) : itemValue}
            onChange={(e) => setItemValue(e.target.value)}
            onClick={e => e.stopPropagation()}
            onKeyDown={e => e.stopPropagation()}
          />
        </Case>
        <Case value={'dataFormat'}>
          <FormControl style={{ marginTop: 10, minWidth: '100px' }}>
            <InputLabel>
              {itemKey}
            </InputLabel>
            <Select
              value={itemValue}
              label="dataFormat"
              size='small'
              onChange={(e: SelectChangeEvent) => setItemValue(e.target.value as string)}
            >
              {Object.entries(dataTypes).map(([_upcasename, name]) => (
                <MenuItem key={name} value={name}>{name}</MenuItem>
              ))}
            </Select>
          </FormControl>
        </Case>
      </SwitchCase>
    </Stack>
  )
}
