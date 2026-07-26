import { DetailsViewItemSingle } from './DetailsViewItemSingle';
import { SelectDataFormat } from './SelectDataFormat';
import { FieldDef } from './fieldDef';
import { DATA_FORMAT_FIELD } from '../../shared/constants';

// The editor for a leaf (non-object) value. `dataFormat` gets a dedicated
// dropdown; everything else uses the generic single-value editor.
export const DetailsViewItemValue = ({
  itemKey,
  itemValue,
  setItemValue,
  def = undefined,
  onCaret = undefined,
}: {
  itemKey: string,
  itemValue: any,
  setItemValue: (newVal: any) => void,
  def?: FieldDef,
  onCaret?: (value: string, caret: number | null, el: HTMLInputElement) => void,
}) => {
  return (
    itemKey === DATA_FORMAT_FIELD
      ? (
        <SelectDataFormat
          value={itemValue}
          setValue={v => setItemValue(v)}
        />
      )
      : (
        <DetailsViewItemSingle
          itemValue={itemValue}
          setItemValue={newVal => setItemValue(newVal)}
          def={def}
          onCaret={onCaret}
        />
      )
  );
};
