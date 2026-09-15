import { MouseEvent, ReactNode } from 'react';
import { DetailsViewItem } from './DetailsViewItem';
import { FieldDef } from './fieldDef';
import { FILES_DIRECTORY_FIELD } from '../../shared/constants';

// The fields of an object row, with the keys sorted as text.
export const DetailsViewItemsInObject = ({
  itemValue,
  setItemValue,
  parentKey,
  def,
  renderBeforeName,
  onRowContextMenu,
  onValueCaret,
  isDescRoot = false,
}: {
  itemValue: { [key: string]: any },
  setItemValue: (newVal: any) => void,
  parentKey: string,
  def?: FieldDef,
  renderBeforeName?: (itemKey: string, parentKey: string | undefined, def?: FieldDef) => ReactNode,
  onRowContextMenu?: (itemKey: string, parentKey: string | undefined, event: MouseEvent<HTMLElement>) => void,
  onValueCaret?: (itemKey: string, parentKey: string | undefined, value: string, caret: number | null, el: HTMLInputElement) => void,
  // True on the document's top-level desc, whose files directory is fixed.
  isDescRoot?: boolean,
}) => {
  return (
    <>
      {Object.entries(itemValue).sort().map(([k, v]) => {
        const isDir = isDescRoot && k === FILES_DIRECTORY_FIELD;

        const changeKey = (newKey: string | undefined) => {
          const item = { ...itemValue };
          delete item[k];
          if (newKey !== undefined) {
            item[newKey] = v;
          }
          setItemValue(item);
        };

        return (
          <DetailsViewItem
            key={k}
            itemKey={k}
            itemValue={v}
            setItemValue={newVal => setItemValue({ ...itemValue, [k]: newVal })}
            setItemKey={isDir ? undefined : changeKey}
            parentKey={parentKey}
            def={def?.children?.[k]}
            renderBeforeName={renderBeforeName}
            onRowContextMenu={onRowContextMenu}
            onValueCaret={onValueCaret}
          />
        );
      })}
    </>
  );
};
