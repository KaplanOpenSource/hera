import { TextField } from '@mui/material';
import { useEffect, useState } from 'react';
import { SelectProperty } from '../../elements/SelectProperty';

enum ItemTypesEnum {
  number = 'number',
  string = 'string',
  null = 'null',
}

const calcItemType = (val: any) => {
  if (val === null) {
    return ItemTypesEnum.null;
  } else if ((typeof val === 'number' || typeof val === 'bigint') && Number.isFinite(val)) {
    return ItemTypesEnum.number;
  } else {
    return ItemTypesEnum.string;
  }
}

export const DetailsViewItemSingle = ({
  itemValue,
  setItemValue,
}: {
  itemValue: any,
  setItemValue: (newVal: any) => void,
}) => {
  const [itemType, setItemType] = useState<ItemTypesEnum>(() => calcItemType(itemValue));

  useEffect(() => {
    setItemType(calcItemType(itemValue));
  }, [itemValue])

  return (
    <>
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
        disabled={itemType === ItemTypesEnum.null}
      />
      <SelectProperty
        label="Type"
        value={itemType}
        setValue={v => setItemType(v as ItemTypesEnum)}
        menuItems={Object.keys(ItemTypesEnum).map((name) => ({ name }))}
      />
    </>
  )
}
