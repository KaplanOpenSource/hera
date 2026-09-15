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
    expect(calcItemType([])).toBe(ItemTypesEnum.array);
    expect(calcItemType([1, 2])).toBe(ItemTypesEnum.array);
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

  it('colors each type chip with a semantic palette color', () => {
    renderTree({ s: 'x', n: 5, o: {}, z: null });
    const chipClass = (label: string) => screen.getByText(label).closest('.MuiChip-root')!.className;
    expect(chipClass('string')).toContain('colorSuccess');
    expect(chipClass('number')).toContain('colorInfo');
    expect(chipClass('null')).toContain('colorWarning');
    // "object" shows on both the parent and child o; at least one carries secondary
    expect(
      screen.getAllByText('object').some(el => el.closest('.MuiChip-root')!.className.includes('colorSecondary'))
    ).toBe(true);
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

describe('list as a field type', () => {
  it('shows the array type and numeric order on a list field', () => {
    render(
      <SimpleTreeView defaultExpandedItems={['config', 'config/a']}>
        <DetailsViewItem itemKey='config' itemValue={{ a: ['x', 'y', 'z'] }} setItemValue={vi.fn()} parentKey={undefined} />
      </SimpleTreeView>
    );
    expect(screen.getByText('array')).toBeDefined();
    const names = screen.getAllByText(/^\d+$/).map(el => el.textContent);
    expect(names).toEqual(['0', '1', '2']);
  });

  it('sorts indices numerically, not as text', () => {
    const many = Array.from({ length: 11 }, (_, i) => 'v' + i);
    render(
      <SimpleTreeView defaultExpandedItems={['config', 'config/a']}>
        <DetailsViewItem itemKey='config' itemValue={{ a: many }} setItemValue={vi.fn()} parentKey={undefined} />
      </SimpleTreeView>
    );
    const names = screen.getAllByText(/^\d+$/).map(el => el.textContent);
    expect(names).toEqual(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']);
  });

  it('keeps a list a list when an element is edited', () => {
    const setItemValue = vi.fn();
    render(
      <SimpleTreeView defaultExpandedItems={['config', 'config/a']}>
        <DetailsViewItem itemKey='config' itemValue={{ a: ['x', 'y'] }} setItemValue={setItemValue} parentKey={undefined} />
      </SimpleTreeView>
    );
    const boxes = screen.getAllByRole('textbox');
    fireEvent.change(boxes[boxes.length - 1], { target: { value: 'changed' } });
    expect(setItemValue).toHaveBeenCalledWith({ a: ['x', 'changed'] });
  });

  it('appends when adding an item to a list', () => {
    const setItemValue = vi.fn();
    render(
      <SimpleTreeView defaultExpandedItems={['config', 'config/a']}>
        <DetailsViewItem itemKey='config' itemValue={{ a: ['x'] }} setItemValue={setItemValue} parentKey={undefined} />
      </SimpleTreeView>
    );
    // The tooltip label sits on a wrapper box, so click the button inside.
    const adds = screen.getAllByLabelText('Add item');
    fireEvent.click(adds[adds.length - 1].querySelector('button')!);
    expect(setItemValue).toHaveBeenCalledWith({ a: ['x', ''] });
  });

  it('removes an element and shifts the rest on delete', () => {
    const setItemValue = vi.fn();
    render(
      <SimpleTreeView defaultExpandedItems={['config', 'config/a']}>
        <DetailsViewItem itemKey='config' itemValue={{ a: ['x', 'y', 'z'] }} setItemValue={setItemValue} parentKey={undefined} />
      </SimpleTreeView>
    );
    fireEvent.click(screen.getByLabelText('Delete 1').querySelector('button')!);
    expect(setItemValue).toHaveBeenCalledWith({ a: ['x', 'z'] });
  });

  it('does not let an index be renamed', () => {
    render(
      <SimpleTreeView defaultExpandedItems={['config', 'config/a']}>
        <DetailsViewItem itemKey='config' itemValue={{ a: ['x'] }} setItemValue={vi.fn()} parentKey={undefined} />
      </SimpleTreeView>
    );
    const before = screen.getAllByRole('textbox').length;
    fireEvent.click(screen.getByText('0'));
    expect(screen.getAllByRole('textbox').length).toBe(before);
  });

  it('converts a list to an object and back, keeping the values', () => {
    const setItemValue = vi.fn();
    render(
      <SimpleTreeView defaultExpandedItems={['config']}>
        <DetailsViewItem itemKey='config' itemValue={{ a: ['x', 'y'] }} setItemValue={setItemValue} parentKey={undefined} />
      </SimpleTreeView>
    );
    pickType('array', 'object');
    expect(setItemValue).toHaveBeenCalledWith({ a: { '0': 'x', '1': 'y' } });
    cleanup();

    const setItemValue2 = vi.fn();
    render(
      <SimpleTreeView defaultExpandedItems={['config']}>
        <DetailsViewItem itemKey='config' itemValue={{ a: { p: 'x', q: 'y' } }} setItemValue={setItemValue2} parentKey={undefined} />
      </SimpleTreeView>
    );
    pickType('object', 'array');
    expect(setItemValue2).toHaveBeenCalledWith({ a: ['x', 'y'] });
  });
});

describe('converting to and from a list', () => {
  it('keeps the order the tree shows when an object becomes a list', () => {
    const setItemValue = vi.fn();
    renderTree({ a: { b: 'B', a: 'A' } }, setItemValue);
    pickType('object', 'array');
    expect(setItemValue).toHaveBeenCalledWith({ a: ['A', 'B'] });
  });

  it('gives zero, not the first element, when a list becomes a number', () => {
    const setItemValue = vi.fn();
    renderTree({ a: ['5', '6'] }, setItemValue);
    pickType('array', 'number');
    expect(setItemValue).toHaveBeenCalledWith({ a: 0 });
  });

  it('gives empty text when a list becomes a string', () => {
    const setItemValue = vi.fn();
    renderTree({ a: ['5', '6'] }, setItemValue);
    pickType('array', 'string');
    expect(setItemValue).toHaveBeenCalledWith({ a: '' });
  });
});

describe('reordering a list by dragging its index', () => {
  it('gives every element a drag handle on its index', () => {
    render(
      <SimpleTreeView defaultExpandedItems={['config', 'config/a']}>
        <DetailsViewItem itemKey='config' itemValue={{ a: ['x', 'y', 'z'] }} setItemValue={vi.fn()} parentKey={undefined} />
      </SimpleTreeView>
    );
    for (const index of [0, 1, 2]) {
      const handle = screen.getByTestId('list-index-' + index);
      expect(handle.getAttribute('role')).toBe('button');
      expect(handle.getAttribute('aria-roledescription')).toBe('sortable');
    }
  });
});
