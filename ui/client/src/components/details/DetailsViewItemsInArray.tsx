import { MouseEvent, ReactNode } from 'react';
import {
  DndContext,
  DragEndEvent,
  KeyboardSensor,
  PointerSensor,
  closestCenter,
  useSensor,
  useSensors,
} from '@dnd-kit/core';
import {
  SortableContext,
  arrayMove,
  sortableKeyboardCoordinates,
  verticalListSortingStrategy,
} from '@dnd-kit/sortable';
import { restrictToVerticalAxis } from '@dnd-kit/modifiers';
import { DetailsViewSortableItem } from './DetailsViewSortableItem';
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
  const sensors = useSensors(
    // A few pixels of movement start a drag, so a plain click still works.
    useSensor(PointerSensor, { activationConstraint: { distance: 4 } }),
    useSensor(KeyboardSensor, { coordinateGetter: sortableKeyboardCoordinates }),
  );

  // Ids are positions, which is all a list element has to identify it by.
  const ids = itemValue.map((_, index) => `${parentKey}#${index}`);

  const onDragEnd = (event: DragEndEvent) => {
    const { active, over } = event;
    if (over && active.id !== over.id) {
      setItemValue(arrayMove(itemValue, ids.indexOf(active.id as string), ids.indexOf(over.id as string)));
    }
  };

  return (
    <DndContext
      sensors={sensors}
      collisionDetection={closestCenter}
      onDragEnd={onDragEnd}
      modifiers={[restrictToVerticalAxis]}
    >
      <SortableContext items={ids} strategy={verticalListSortingStrategy}>
        {itemValue.map((value, index) => (
          <DetailsViewSortableItem
            key={ids[index]}
            id={ids[index]}
            index={index}
            itemValue={value}
            setItemValue={newVal => setItemValue(itemValue.map((v, i) => (i === index ? newVal : v)))}
            // An index can't be renamed, so only deletion goes through here.
            setItemKey={newKey => {
              if (newKey === undefined) {
                setItemValue(itemValue.filter((_, i) => i !== index));
              }
            }}
            parentKey={parentKey}
            def={def}
            renderBeforeName={renderBeforeName}
            onRowContextMenu={onRowContextMenu}
            onValueCaret={onValueCaret}
          />
        ))}
      </SortableContext>
    </DndContext>
  );
};
