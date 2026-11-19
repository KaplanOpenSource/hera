import { describe, it, expect, vi, afterEach } from 'vitest';
import { render, screen, fireEvent, cleanup } from '@testing-library/react';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { DetailsViewItem } from '../src/components/details/DetailsViewItem';
import { calcItemType, coerceToType, ItemTypesEnum } from '../src/components/details/ItemTypeSelector';

const renderTree = (itemValue: any, setItemValue = vi.fn()) => {
  render(
    <SimpleTreeView defaultExpandedItems={['config']}>
      <DetailsViewItem itemKey='config' itemValue={itemValue} setItemValue={setItemValue} parentKey={undefined} />
    </SimpleTreeView>
  );
  return setItemValue;
};

// Open the child field's type chip menu and pick a type.
const pickType = (currentLabel: string, target: string) => {
  const chips = screen.getAllByText(currentLabel);
  fireEvent.click(chips[chips.length - 1]);
  fireEvent.click(screen.getByRole('menuitem', { name: target }));
};

afterEach(() => cleanup());

describe('calcItemType for booleans', () => {
  it('classifies true and false as boolean, not string', () => {
    expect(calcItemType(true)).toBe(ItemTypesEnum.boolean);
    expect(calcItemType(false)).toBe(ItemTypesEnum.boolean);
    expect(calcItemType('true')).toBe(ItemTypesEnum.string);
  });
});

describe('coerceToType for booleans', () => {
  it('turns values into a boolean', () => {
    expect(coerceToType(ItemTypesEnum.boolean, 'true')).toBe(true);
    expect(coerceToType(ItemTypesEnum.boolean, 'TRUE')).toBe(true);
    expect(coerceToType(ItemTypesEnum.boolean, 'yes')).toBe(false);
    expect(coerceToType(ItemTypesEnum.boolean, '')).toBe(false);
    expect(coerceToType(ItemTypesEnum.boolean, 5)).toBe(true);
    expect(coerceToType(ItemTypesEnum.boolean, 0)).toBe(false);
    expect(coerceToType(ItemTypesEnum.boolean, null)).toBe(false);
    expect(coerceToType(ItemTypesEnum.boolean, {})).toBe(false);
    expect(coerceToType(ItemTypesEnum.boolean, [])).toBe(false);
    expect(coerceToType(ItemTypesEnum.boolean, true)).toBe(true);
  });

  it('turns a boolean into text and a number', () => {
    expect(coerceToType(ItemTypesEnum.string, true)).toBe('true');
    expect(coerceToType(ItemTypesEnum.string, false)).toBe('false');
    expect(coerceToType(ItemTypesEnum.number, true)).toBe(1);
    expect(coerceToType(ItemTypesEnum.number, false)).toBe(0);
  });
});

describe('boolean field editor', () => {
  it('shows a TRUE switch, not a text box', () => {
    renderTree({ a: true });
    expect(screen.getByText('boolean')).toBeDefined();
    expect(screen.getByText('TRUE')).toBeDefined();
    expect((screen.getByRole('switch') as HTMLInputElement).checked).toBe(true);
    expect(screen.queryByDisplayValue('true')).toBeNull();
  });

  it('shows FALSE for a false value', () => {
    renderTree({ a: false });
    expect(screen.getByText('FALSE')).toBeDefined();
    expect((screen.getByRole('switch') as HTMLInputElement).checked).toBe(false);
  });

  it('saves a real boolean when toggled', () => {
    const setItemValue = renderTree({ a: true });
    fireEvent.click(screen.getByRole('switch'));
    expect(setItemValue).toHaveBeenCalledWith({ a: false });
  });

  it('converts a "true" string to a boolean from the type menu', () => {
    const setItemValue = renderTree({ a: 'true' });
    pickType('string', 'boolean');
    expect(setItemValue).toHaveBeenCalledWith({ a: true });
  });

  it('colors the boolean chip', () => {
    renderTree({ a: true });
    expect(screen.getByText('boolean').closest('.MuiChip-root')!.className).toContain('colorError');
  });
});
