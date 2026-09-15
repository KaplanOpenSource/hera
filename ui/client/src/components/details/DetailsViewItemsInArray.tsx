import { MouseEvent, ReactNode } from 'react';
import { DetailsViewItem } from './DetailsViewItem';
import { FieldDef } from './fieldDef';

// The elements of a list row, kept in their own order and always still a list.
export const DetailsViewItemsInArray = ({
  itemValue,
  setItemValue,
  parentKey,
  def,
  renderBeforeName,
  onRowContextMenu,
  onValueCaret,
}: {
  itemValue: any[],
  setItemValue: (newVal: any) => void,
  parentKey: string,
  def?: FieldDef,
  renderBeforeName?: (itemKey: string, parentKey: string | undefined, def?: FieldDef) => ReactNode,
  onRowContextMenu?: (itemKey: string, parentKey: string | undefined, event: MouseEvent<HTMLElement>) => void,
  onValueCaret?: (itemKey: string, parentKey: string | undefined, value: string, caret: number | null, el: HTMLInputElement) => void,
}) => {
  return (
    <>
      {itemValue.map((value, index) => {
        const setElement = (newVal: any) => {
          setItemValue(itemValue.map((v, i) => (i === index ? newVal : v)));
        };

        // An index can't be renamed, so only deletion goes through here.
        const deleteElement = (newKey: string | undefined) => {
          if (newKey === undefined) {
            setItemValue(itemValue.filter((_, i) => i !== index));
          }
        };

        return (
          <DetailsViewItem
            key={index}
            itemKey={String(index)}
            itemValue={value}
            setItemValue={setElement}
            setItemKey={deleteElement}
            isListIndex
            parentKey={parentKey}
            def={def?.children?.[String(index)]}
            renderBeforeName={renderBeforeName}
            onRowContextMenu={onRowContextMenu}
            onValueCaret={onValueCaret}
          />
        );
      })}
    </>
  );
};
