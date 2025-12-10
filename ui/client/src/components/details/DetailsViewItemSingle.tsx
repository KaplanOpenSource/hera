import { FormControl, InputLabel, MenuItem, Select, Stack, TextField, Typography } from '@mui/material';
import { Case, SwitchCase } from '../../elements/SwitchCase';
import { useServerConstants } from '../../stores/useServerConstants';

export const DetailsViewItemSingle = ({
  itemKey,
  itemValue,
}: {
  itemKey: string,
  itemValue: any,
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
            value={
              JSON.stringify(itemValue)
            }
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
            // onChange={(event: SelectChangeEvent) => setAge(event.target.value as string)}
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
