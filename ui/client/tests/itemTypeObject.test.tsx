import { describe, it, expect, vi, afterEach } from 'vitest';
import { render, screen, fireEvent, cleanup } from '@testing-library/react';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { DetailsViewItem } from '../src/components/details/DetailsViewItem';
import { calcItemType, ItemTypesEnum } from '../src/components/details/ItemTypeSelector';

const renderTree = (itemValue: any, setItemValue = vi.fn()) => {
  render(
    <SimpleTreeView defaultExpandedItems={['config']}>
      <DetailsViewItem itemKey='config' itemValue={itemValue} setItemValue={setItemValue} parentKey={undefined} />
    </SimpleTreeView>
  );
  return setItemValue;
};

// Open the child field's type chip menu and pick a type. The parent row also
// has a chip, so target the last chip with the given label (the child's).
const pickType = (currentLabel: string, target: string) => {
  const chips = screen.getAllByText(currentLabel);
  fireEvent.click(chips[chips.length - 1]);
  fireEvent.click(screen.getByRole('menuitem', { name: target }));
};

afterEach(() => cleanup());

describe('calcItemType', () => {
  it('classifies objects, arrays, null, numbers and strings', () => {
    expect(calcItemType({})).toBe(ItemTypesEnum.object);
    expect(calcItemType({ a: 1 })).toBe(ItemTypesEnum.object);
    expect(calcItemType([])).toBe(ItemTypesEnum.object);
    expect(calcItemType(null)).toBe(ItemTypesEnum.null);
    expect(calcItemType(42)).toBe(ItemTypesEnum.number);
    expect(calcItemType('hi')).toBe(ItemTypesEnum.string);
  });
});

describe('object as a field type', () => {
  it('turns a scalar field into an empty substructure when set to object', () => {
    const setItemValue = renderTree({ a: 'x' });
    pickType('string', 'object');
    // child a's setter rewrites the parent object with a now {}
    expect(setItemValue).toHaveBeenCalledWith({ a: {} });
  });

  it('reverts an object field back to text when set to string', () => {
    const setItemValue = renderTree({ a: { nested: '1' } });
    pickType('object', 'string');
    expect(setItemValue).toHaveBeenCalledWith({ a: '' });
  });

  it('shows the object type on an object-valued field', () => {
    renderTree({ a: { nested: '1' } });
    // both the parent (config) and child (a) are objects
    expect(screen.getAllByText('object').length).toBe(2);
  });
});
