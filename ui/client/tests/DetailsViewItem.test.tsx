import { describe, it, expect, vi, afterEach } from 'vitest';
import { render, screen, cleanup } from '@testing-library/react';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { DetailsViewItem } from '../src/components/details/DetailsViewItem';

const renderItem = ({
  itemKey,
  itemValue,
  parentKey,
}: {
  itemKey: string;
  itemValue: any;
  parentKey?: string;
}) => {
  return render(
    <SimpleTreeView defaultExpandedItems={[parentKey ? `${parentKey}/${itemKey}` : itemKey]}>
      <DetailsViewItem
        itemKey={itemKey}
        itemValue={itemValue}
        setItemValue={vi.fn()}
        parentKey={parentKey}
      />
    </SimpleTreeView>
  );
};

describe('DetailsViewItem', () => {
  afterEach(() => {
    cleanup();
  });

  it('renders a leaf item with its key name', () => {
    renderItem({ itemKey: 'name', itemValue: 'hello' });
    expect(screen.getByText('name')).toBeDefined();
  });

  it('renders a leaf item with its string value', () => {
    renderItem({ itemKey: 'name', itemValue: 'hello' });
    expect(screen.getByDisplayValue('hello')).toBeDefined();
  });

  it('renders an object item with its key name', () => {
    renderItem({ itemKey: 'config', itemValue: { a: '1' } });
    expect(screen.getByText('config')).toBeDefined();
  });

  it('renders sub-items of an object', () => {
    renderItem({ itemKey: 'config', itemValue: { myField: 'myValue' } });
    expect(screen.getByText('myField')).toBeDefined();
  });

  it('shows "(empty)" for an empty object', () => {
    renderItem({ itemKey: 'emptyObj', itemValue: {} });
    expect(screen.getByText('(empty)')).toBeDefined();
  });

  it('does not show "(empty)" for a non-empty object', () => {
    renderItem({ itemKey: 'config', itemValue: { a: '1' } });
    expect(screen.queryByText('(empty)')).toBeNull();
  });

  it('does not show "(empty)" for a leaf value', () => {
    renderItem({ itemKey: 'name', itemValue: 'hello' });
    expect(screen.queryByText('(empty)')).toBeNull();
  });
});
