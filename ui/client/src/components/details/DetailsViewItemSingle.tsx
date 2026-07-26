import { TextField } from '@mui/material';
import { useEffect, useState } from 'react';
import { FieldDef } from './fieldDef';
import { ItemTypesEnum, calcItemType } from './ItemTypeSelector';

export const DetailsViewItemSingle = ({
  itemValue,
  setItemValue,
  def = undefined,
  onCaret = undefined,
}: {
  itemValue: any,
  setItemValue: (newVal: any) => void,
  def?: FieldDef,
  // Optional: reports this field's value and caret position after any edit or
  // caret move (typing, arrow keys, clicks) so a caller can drive inline
  // autocomplete. Not provided elsewhere in the details view.
  onCaret?: (value: string, caret: number | null, el: HTMLInputElement) => void,
}) => {
  const [itemType, setItemType] = useState<ItemTypesEnum>(() => calcItemType(itemValue));
  const missing = !!def?.required && (itemValue === undefined || itemValue === null || itemValue === '');

  useEffect(() => {
    setItemType(calcItemType(itemValue));
  }, [itemValue])

  const reportCaret = (el: HTMLInputElement) => {
    onCaret?.(el.value, el.selectionStart, el);
  };

  return (
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
        reportCaret(e.target as HTMLInputElement);
      }}
      onClick={e => { e.stopPropagation(); reportCaret(e.currentTarget.querySelector('input') ?? e.target as HTMLInputElement); }}
      onKeyUp={e => reportCaret(e.target as HTMLInputElement)}
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
  )
}
