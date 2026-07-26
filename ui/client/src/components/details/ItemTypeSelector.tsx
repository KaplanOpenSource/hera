import { useState, MouseEvent } from 'react';
import { Chip, Menu, MenuItem } from '@mui/material';

export enum ItemTypesEnum {
  number = 'number',
  string = 'string',
  null = 'null',
}

// The type a raw value currently holds — drives which editor is shown.
export const calcItemType = (val: any) => {
  if (val === null) {
    return ItemTypesEnum.null;
  } else if ((typeof val === 'number' || typeof val === 'bigint') && Number.isFinite(val)) {
    return ItemTypesEnum.number;
  } else {
    return ItemTypesEnum.string;
  }
};

// A small chip showing the value's type. Clicking it opens a menu to pick a
// type, which both records the new type and coerces the current value to match
// (null / number / string).
export const ItemTypeSelector = ({
  itemType,
  setItemType,
  itemValue,
  setItemValue,
}: {
  itemType: ItemTypesEnum,
  setItemType: (t: ItemTypesEnum) => void,
  itemValue: any,
  setItemValue: (newVal: any) => void,
}) => {
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);

  const chooseType = (v: ItemTypesEnum) => {
    setItemType(v);
    if (v === ItemTypesEnum.null) {
      setItemValue(null);
    } else if (v === ItemTypesEnum.number) {
      const num = parseFloat(itemValue);
      setItemValue(Number.isFinite(num) ? num : 0);
    } else {
      setItemValue(itemValue + '');
    }
    setAnchorEl(null);
  };

  return (
    <>
      <Chip
        size="small"
        variant="outlined"
        label={itemType}
        onClick={(e: MouseEvent<HTMLElement>) => {
          e.stopPropagation();
          setAnchorEl(e.currentTarget);
        }}
      />
      <Menu
        anchorEl={anchorEl}
        open={!!anchorEl}
        onClose={() => setAnchorEl(null)}
      >
        {Object.values(ItemTypesEnum).map((t) => (
          <MenuItem
            key={t}
            selected={t === itemType}
            onClick={() => chooseType(t)}
          >
            {t}
          </MenuItem>
        ))}
      </Menu>
    </>
  );
};
