import { useState, MouseEvent } from 'react';
import { Chip, Menu, MenuItem } from '@mui/material';

export enum ItemTypesEnum {
  number = 'number',
  string = 'string',
  null = 'null',
  object = 'object',
}

// Each type gets a semantic MUI palette color so it adapts to light/dark mode.
const typeColor: { [key in ItemTypesEnum]: 'success' | 'info' | 'secondary' | 'warning' } = {
  [ItemTypesEnum.string]: 'success',
  [ItemTypesEnum.number]: 'info',
  [ItemTypesEnum.object]: 'secondary',
  [ItemTypesEnum.null]: 'warning',
};

// The type a raw value currently holds — drives which editor is shown.
export const calcItemType = (val: any) => {
  if (val === null) {
    return ItemTypesEnum.null;
  } else if (typeof val === 'object') {
    return ItemTypesEnum.object;
  } else if ((typeof val === 'number' || typeof val === 'bigint') && Number.isFinite(val)) {
    return ItemTypesEnum.number;
  } else {
    return ItemTypesEnum.string;
  }
};

// The value to store when a field is switched to the given type. Objects become
// an empty substructure; scalars coerce the current value where it makes sense.
const coerceToType = (t: ItemTypesEnum, current: any) => {
  switch (t) {
    case ItemTypesEnum.object:
      return {};
    case ItemTypesEnum.null:
      return null;
    case ItemTypesEnum.number: {
      const num = parseFloat(current);
      return Number.isFinite(num) ? num : 0;
    }
    default: {
      // string: an object/null has no sensible text, so start empty.
      const currentType = calcItemType(current);
      return (currentType === ItemTypesEnum.null || currentType === ItemTypesEnum.object) ? '' : String(current);
    }
  }
};

// A small chip showing the value's type. Clicking it opens a menu to pick a
// type, which coerces the current value to match (including object = substructure).
export const ItemTypeSelector = ({
  itemValue,
  setItemValue,
}: {
  itemValue: any,
  setItemValue: (newVal: any) => void,
}) => {
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);
  const itemType = calcItemType(itemValue);

  const chooseType = (t: ItemTypesEnum) => {
    setItemValue(coerceToType(t, itemValue));
    setAnchorEl(null);
  };

  return (
    <>
      <Chip
        size="small"
        variant="outlined"
        color={typeColor[itemType]}
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
            sx={{ color: theme => theme.palette[typeColor[t]].main }}
          >
            {t}
          </MenuItem>
        ))}
      </Menu>
    </>
  );
};
