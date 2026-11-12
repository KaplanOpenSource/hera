import { useState, MouseEvent } from 'react';
import { Chip, Menu, MenuItem } from '@mui/material';

export enum ItemTypesEnum {
  number = 'number',
  string = 'string',
  null = 'null',
  object = 'object',
  array = 'array',
}

// Each type gets a semantic MUI palette color so it adapts to light/dark mode.
const typeColor: { [key in ItemTypesEnum]: 'success' | 'info' | 'secondary' | 'warning' | 'primary' } = {
  [ItemTypesEnum.string]: 'success',
  [ItemTypesEnum.number]: 'info',
  [ItemTypesEnum.object]: 'secondary',
  [ItemTypesEnum.array]: 'primary',
  [ItemTypesEnum.null]: 'warning',
};

// The type a raw value currently holds — drives which editor is shown.
export const calcItemType = (val: any) => {
  if (val === null) {
    return ItemTypesEnum.null;
  } else if (Array.isArray(val)) {
    return ItemTypesEnum.array;
  } else if (typeof val === 'object') {
    return ItemTypesEnum.object;
  } else if ((typeof val === 'number' || typeof val === 'bigint') && Number.isFinite(val)) {
    return ItemTypesEnum.number;
  } else {
    return ItemTypesEnum.string;
  }
};

// Types that hold no single value to carry over into a scalar.
const noValueTypes = [ItemTypesEnum.null, ItemTypesEnum.object, ItemTypesEnum.array];

// The value to store when a field is switched to the given type.
const coerceToType = (t: ItemTypesEnum, current: any) => {
  switch (t) {
    case ItemTypesEnum.object:
      if (Array.isArray(current)) {
        return Object.fromEntries(current.map((v, i) => [String(i), v]));
      }
      return {};
    case ItemTypesEnum.array:
      if (Array.isArray(current)) {
        return current;
      }
      if (calcItemType(current) === ItemTypesEnum.object) {
        // Sorted like the tree shows them, so the order matches what was seen.
        return Object.entries(current).sort().map(([, v]) => v);
      }
      return [];
    case ItemTypesEnum.null:
      return null;
    case ItemTypesEnum.number: {
      // parseFloat on a list would wrongly pick up its first element.
      if (noValueTypes.includes(calcItemType(current))) {
        return 0;
      }
      const num = parseFloat(current);
      return Number.isFinite(num) ? num : 0;
    }
    default: {
      // string: an object or a list has no sensible text, so start empty.
      return noValueTypes.includes(calcItemType(current)) ? '' : String(current);
    }
  }
};

// A chip showing the value's type. Click it to pick another type.
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
