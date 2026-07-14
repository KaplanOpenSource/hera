import { TextField } from '@mui/material';
import { useEffect, useState } from 'react';
import { FieldDef } from './fieldDef';
import { ItemTypeSelector, ItemTypesEnum, calcItemType } from './ItemTypeSelector';

export const DetailsViewItemSingle = ({
  itemValue,
  setItemValue,
  def = undefined,
}: {
  itemValue: any,
  setItemValue: (newVal: any) => void,
  def?: FieldDef,
}) => {
  const [itemType, setItemType] = useState<ItemTypesEnum>(() => calcItemType(itemValue));
  const missing = !!def?.required && (itemValue === undefined || itemValue === null || itemValue === '');

  useEffect(() => {
    setItemType(calcItemType(itemValue));
  }, [itemValue])

  return (
    <>
      <TextField
        size='small'
        value={['object', 'function'].includes(typeof itemValue) ? JSON.stringify(itemValue) : itemValue}
        onChange={(e) => {
          if (itemType === ItemTypesEnum.number) {
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
        error={missing}
        helperText={missing ? 'required' : undefined}
        // Float the "required" text just below the field: absolute so it adds no
        // height (required and normal rows stay the same height), zIndex so it is
        // not painted over by the next row.
        sx={{ position: 'relative' }}
        slotProps={{
          formHelperText: {
            sx: { position: 'absolute', top: '100%', left: 0, m: 0, lineHeight: 1, zIndex: 1 },
          },
        }}
      />
      <ItemTypeSelector
        itemType={itemType}
        setItemType={setItemType}
        itemValue={itemValue}
        setItemValue={setItemValue}
      />
    </>
  )
}
