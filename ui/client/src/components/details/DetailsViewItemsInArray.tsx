import { MouseEvent, ReactNode } from 'react';
import { DetailsViewItem } from './DetailsViewItem';
import { DetailsViewListIndex } from './DetailsViewListIndex';
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
  // Drag an index onto another one to move that element there.
  const moveElement = (from: number, to: number) => {
    if (from === to) {
      return;
    }
    const next = [...itemValue];
    const [moved] = next.splice(from, 1);
    next.splice(to, 0, moved);
    setItemValue(next);
  };

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
            nameView={<DetailsViewListIndex index={index} onMove={moveElement} />}
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
