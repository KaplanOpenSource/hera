import { SelectProperty } from '../../elements/SelectProperty';

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

// Dropdown that picks a value's type. Changing the type both records the new
// type and coerces the current value to match (null / number / string).
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
  return (
    <SelectProperty
      label="Type"
      value={itemType}
      setValue={v => {
        setItemType(v as ItemTypesEnum);
        if (v === ItemTypesEnum.null) {
          setItemValue(null);
        } else if (v === ItemTypesEnum.number) {
          const num = parseFloat(itemValue);
          setItemValue(Number.isFinite(num) ? num : 0);
        } else {
          setItemValue(itemValue + '');
        }
      }}
      menuItems={Object.keys(ItemTypesEnum).map((name) => ({ name }))}
    />
  );
};
