import { describe, it, expect, vi, afterEach } from 'vitest';
import { useState } from 'react';
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

  it('shows no value editor for a null field (just the type chip)', () => {
    render(
      <SimpleTreeView defaultExpandedItems={['x']}>
        <DetailsViewItem itemKey='x' itemValue={null} setItemValue={vi.fn()} parentKey={undefined} />
      </SimpleTreeView>
    );
    expect(screen.getByText('null')).toBeDefined();      // the type chip
    expect(screen.queryByRole('textbox')).toBeNull();    // no disabled "null" field
  });

  it('expands a field when it becomes an object', () => {
    // Stateful harness so the value actually turns into {} and the row can open.
    const Harness = () => {
      const [val, setVal] = useState<any>({ a: 'x' });
      return (
        <SimpleTreeView defaultExpandedItems={['config']}>
          <DetailsViewItem itemKey='config' itemValue={val} setItemValue={setVal} parentKey={undefined} />
        </SimpleTreeView>
      );
    };
    render(<Harness />);
    // a is a leaf, so its "(empty)" substructure label is not shown yet
    expect(screen.queryByText('(empty)')).toBeNull();
    pickType('string', 'object');
    // a is now an empty object AND auto-expanded, so its "(empty)" label shows
    expect(screen.getByText('(empty)')).toBeDefined();
  });
});
