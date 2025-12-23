import { FormControl, InputLabel, MenuItem, Select, SelectChangeEvent, Stack, TextField, Typography } from '@mui/material';
import { Case, SwitchCase } from '../../elements/SwitchCase';
import { useServerConstants } from '../../stores/useServerConstants';
import { useState } from 'react';

enum ItemTypesEnum {
  number = 'number',
  string = 'string',
}

export const DetailsViewItemSingle = ({
  itemKey,
  itemValue,
  setItemValue,
}: {
  itemKey: string,
  itemValue: any,
  setItemValue: (newVal: any) => void,
}) => {
  const { dataTypes } = useServerConstants();
  const [itemType, setItemType] = useState<ItemTypesEnum>((Number.isFinite(itemValue) ? 'number' : 'string') as ItemTypesEnum);

  return (
    <Stack direction='row' spacing={1} justifyItems={'center'} alignItems={'center'} style={{ marginTop: 8 }}>
      <Typography>
        {itemKey}
      </Typography>
      <SwitchCase test={itemKey}>
        <Case isDefault>
          <TextField
            size='small'
            value={['object', 'function'].includes(typeof itemValue) ? JSON.stringify(itemValue) : itemValue}
            onChange={(e) => {
              if (itemType === 'number') {
                setItemValue(parseFloat(e.target.value))
              } else {
                setItemValue(e.target.value)
              }
            }}
            onClick={e => e.stopPropagation()}
            onKeyDown={e => e.stopPropagation()}
          />
          <FormControl style={{ minWidth: '100px' }}>
            <InputLabel>
              Type
            </InputLabel>
            <Select
              value={itemType}
              label="Type"
              size='small'
              onChange={e => setItemType(e.target.value)}
            >
              {Object.keys(ItemTypesEnum).filter(x => isNaN(x as any)).map(name => (
                <MenuItem key={name} value={name}>{name}</MenuItem>
              ))}
            </Select>
          </FormControl>
        </Case>
        <Case value={'dataFormat'}>
          <FormControl style={{ minWidth: '100px' }}>
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
