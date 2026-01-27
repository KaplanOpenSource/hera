import { FormControl, InputLabel, MenuItem, Select, TextField } from '@mui/material';
import { useState } from 'react';
import { Case, SwitchCase } from '../../elements/SwitchCase';
import { SelectDataFormat } from './SelectDataFormat';
import { SelectProperty } from '../../elements/SelectProperty';

enum ItemTypesEnum {
  number = 'number',
  string = 'string',
}

const isnum = (v: any) => {
  return ['number', 'bigint'].includes(typeof v) && Number.isFinite(v);
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
  const [itemType, setItemType] = useState<ItemTypesEnum>((isnum(itemValue) ? 'number' : 'string') as ItemTypesEnum);

  return (
    <SwitchCase test={itemKey}>
      <Case isDefault>
        <TextField
          size='small'
          value={['object', 'function'].includes(typeof itemValue) ? JSON.stringify(itemValue) : itemValue}
          onChange={(e) => {
            if (itemType === 'number') {
              const num = parseFloat(e.target.value);
              if (Number.isFinite(num)) {
                setItemValue(num)
              }
            } else {
              setItemValue(e.target.value)
            }
          }}
          onClick={e => e.stopPropagation()}
          onKeyDown={e => e.stopPropagation()}
          fullWidth
        />
        <SelectProperty
          label="Type"
          value={itemType}
          setValue={v => setItemType(v as ItemTypesEnum)}
          menuItems={Object.keys(ItemTypesEnum).map((name) => ({ name }))}
        />
      </Case>
      <Case value={'dataFormat'}>
        <SelectDataFormat
          value={itemValue}
          setValue={v => setItemValue(v)}
        />
      </Case>
    </SwitchCase>
  )
}
