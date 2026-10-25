import { describe, it, expect, vi, afterEach } from 'vitest';
import { render, screen, fireEvent, cleanup, act } from '@testing-library/react';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import { DetailsViewItem } from '../src/components/details/DetailsViewItem';

afterEach(() => cleanup());

// jsdom gives every element a zero-size rect, which dnd-kit needs to detect
// collisions, so give each list row a real one based on its index.
const stubRowSizes = () => {
  const rowHeight = 40;
  const rect = (top: number, height: number) => {
    return { x: 0, y: top, top, left: 0, bottom: top + height, right: 200, width: 200, height, toJSON: () => {} } as DOMRect;
  };
  Element.prototype.getBoundingClientRect = function () {
    const self = this as HTMLElement;
    const handles = self.querySelectorAll?.('[data-testid^="list-index-"]') ?? [];
    // A container spans all the rows it holds; a row spans just its own.
    if (handles.length > 1) {
      return rect(0, handles.length * rowHeight);
    }
    const testid = self.dataset?.testid ?? (handles[0] as HTMLElement | undefined)?.dataset?.testid;
    const index = testid ? Number(testid.replace('list-index-', '')) : 0;
    return rect(index * rowHeight, rowHeight);
  };
};

const renderList = (value: any[]) => {
  stubRowSizes();
  const setItemValue = vi.fn();
  render(
    <SimpleTreeView defaultExpandedItems={['config', 'config/a']}>
      <DetailsViewItem itemKey='config' itemValue={{ a: value }} setItemValue={setItemValue} parentKey={undefined} />
    </SimpleTreeView>
  );
  return setItemValue;
};

// dnd-kit measures between steps, so let a tick pass after each one.
const tick = async () => {
  await act(async () => { await new Promise(r => setTimeout(r, 0)); });
};

// Drive dnd-kit's keyboard sensor: space starts the drag, arrows move, space drops.
const dragDown = async (fromIndex: number, steps: number) => {
  const handle = screen.getByTestId('list-index-' + fromIndex);
  handle.focus();
  fireEvent.keyDown(handle, { code: 'Space', key: ' ' });
  await tick();
  for (let i = 0; i < steps; i++) {
    fireEvent.keyDown(handle, { code: 'ArrowDown', key: 'ArrowDown' });
    await tick();
  }
  fireEvent.keyDown(handle, { code: 'Space', key: ' ' });
  await tick();
};

// Drive dnd-kit's pointer sensor, which is what a mouse drag uses.
const dragWithMouse = async (fromIndex: number, toIndex: number) => {
  const handle = screen.getByTestId('list-index-' + fromIndex);
  fireEvent.pointerDown(handle, { button: 0, isPrimary: true, pointerId: 1, pointerType: 'mouse', clientX: 10, clientY: fromIndex * 40 + 20 });
  await tick();
  fireEvent.pointerMove(document, { isPrimary: true, pointerId: 1, clientX: 10, clientY: fromIndex * 40 + 40 });
  await tick();
  fireEvent.pointerMove(document, { isPrimary: true, pointerId: 1, clientX: 10, clientY: toIndex * 40 + 20 });
  await tick();
  fireEvent.pointerUp(document, { isPrimary: true, pointerId: 1, clientX: 10, clientY: toIndex * 40 + 20 });
  await tick();
};

describe('reordering a list by dragging its index', () => {
  it('moves an element down with the mouse', async () => {
    const setItemValue = renderList(['x', 'y', 'z']);
    await dragWithMouse(0, 1);
    expect(setItemValue).toHaveBeenCalledWith({ a: ['y', 'x', 'z'] });
  });

  it('starts a drag when space is pressed on the handle', async () => {
    renderList(['x', 'y', 'z']);
    const handle = screen.getByTestId('list-index-0');
    handle.focus();
    expect(document.activeElement).toBe(handle);
    fireEvent.keyDown(handle, { code: 'Space', key: ' ' });
    await act(async () => { await new Promise(r => setTimeout(r, 0)); });
    expect(handle.getAttribute('aria-pressed')).toBe('true');
  });

  it('moves an element down', async () => {
    const setItemValue = renderList(['x', 'y', 'z']);
    await dragDown(0, 1);
    expect(setItemValue).toHaveBeenCalledWith({ a: ['y', 'x', 'z'] });
  });
});
