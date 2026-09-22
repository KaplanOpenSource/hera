import { MouseEvent, ReactNode } from 'react';
import { useSortable } from '@dnd-kit/sortable';
import { CSS } from '@dnd-kit/utilities';
import { DetailsViewItem } from './DetailsViewItem';
import { DetailsViewListIndex } from './DetailsViewListIndex';
import { FieldDef } from './fieldDef';

// One list element, which slides aside while another one is dragged over it.
export const DetailsViewSortableItem = ({
  id,
  index,
  itemValue,
  setItemValue,
  setItemKey,
  parentKey,
  def,
  nameForView,
  renderBeforeName,
  onRowContextMenu,
  onValueCaret,
}: {
  id: string,
  index: number,
  itemValue: any,
  setItemValue: (newVal: any) => void,
  setItemKey: (newKey: string | undefined) => void,
  parentKey: string,
  def?: FieldDef,
  nameForView?: (itemKey: string, parentKey: string | undefined) => ReactNode,
  renderBeforeName?: (itemKey: string, parentKey: string | undefined, def?: FieldDef) => ReactNode,
  onRowContextMenu?: (itemKey: string, parentKey: string | undefined, event: MouseEvent<HTMLElement>) => void,
  onValueCaret?: (itemKey: string, parentKey: string | undefined, value: string, caret: number | null, el: HTMLInputElement) => void,
}) => {
  // Ids are positions, so after a drop the order of ids never changes. Without
  // this dnd-kit would animate the dropped row back to where it started.
  const { attributes, listeners, setNodeRef, transform, transition, isDragging } = useSortable({
    id,
    animateLayoutChanges: () => false,
  });

  const style = {
    transform: CSS.Transform.toString(transform),
    transition,
    opacity: 1,
  };
  if (isDragging) {
    style.opacity = 0.4;
  }

  return (
    <DetailsViewItem
      rootRef={setNodeRef}
      rootStyle={style}
      itemKey={String(index)}
      itemValue={itemValue}
      setItemValue={setItemValue}
      setItemKey={setItemKey}
      // A list element is named by its position, which the user can't rename.
      // Rows under it keep the name function passed from above.
      nameForView={(k, pk) => {
        if (pk === parentKey && k === String(index)) {
          return <DetailsViewListIndex index={index} attributes={attributes} listeners={listeners} />;
        }
        return nameForView?.(k, pk);
      }}
      allowRename={false}
      parentKey={parentKey}
      def={def?.children?.[String(index)]}
      renderBeforeName={renderBeforeName}
      onRowContextMenu={onRowContextMenu}
      onValueCaret={onValueCaret}
    />
  );
};
