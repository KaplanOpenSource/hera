import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, fireEvent, render, screen } from '@testing-library/react';
import { WorkflowContextMenu, WorkflowContextMenuKind } from '../src/components/workflow/WorkflowContextMenu';

afterEach(() => cleanup());

const referenceOptions = [
  { node: 'src', outputs: ['out1', 'out2'] },
  { node: 'other', outputs: ['bee'] },
];

describe('WorkflowContextMenu', () => {
  it('renders nothing when there is no target', () => {
    render(<WorkflowContextMenu menu={null} referenceOptions={[]} onClose={vi.fn()} onDeleteNode={vi.fn()} onDeleteField={vi.fn()} onRemoveRequire={vi.fn()} onReferenceOutput={vi.fn()} />);
    expect(screen.queryByText(/Delete node/)).toBeNull();
    expect(screen.queryByText(/Remove requirement/)).toBeNull();
  });

  it('references another node output via two in-menu steps', async () => {
    const onReferenceOutput = vi.fn();
    const onClose = vi.fn();
    render(
      <WorkflowContextMenu
        menu={{ kind: WorkflowContextMenuKind.Field, node: 'alpha', param: 'p', x: 10, y: 20 }}
        referenceOptions={referenceOptions}
        onClose={onClose}
        onDeleteNode={vi.fn()}
        onDeleteField={vi.fn()}
        onRemoveRequire={vi.fn()}
        onReferenceOutput={onReferenceOutput}
      />,
    );
    // Clicking the action switches the menu into step 1 (pick a node) without
    // closing or calling the handler yet.
    fireEvent.click(screen.getByText(/Reference another node's output/));
    expect(onReferenceOutput).not.toHaveBeenCalled();
    expect(onClose).not.toHaveBeenCalled();

    // Step 1: pick the source node (only one autocomplete shows per step).
    fireEvent.change(screen.getByRole('combobox'), { target: { value: 'src' } });
    fireEvent.click(await screen.findByRole('option', { name: 'src' }));

    // Step 2: pick one of that node's outputs → inserts the reference and closes.
    fireEvent.change(screen.getByRole('combobox'), { target: { value: 'out2' } });
    fireEvent.click(await screen.findByRole('option', { name: 'out2' }));
    expect(onReferenceOutput).toHaveBeenCalledWith('alpha', 'p', 'src', 'out2');
    expect(onClose).toHaveBeenCalled();
  });

  it('deletes a node and closes when the node item is clicked', () => {
    const onDeleteNode = vi.fn();
    const onClose = vi.fn();
    render(
      <WorkflowContextMenu
        menu={{ kind: WorkflowContextMenuKind.Node, name: 'alpha', x: 10, y: 20 }}
        referenceOptions={[]}
        onClose={onClose}
        onDeleteNode={onDeleteNode}
        onDeleteField={vi.fn()}
        onRemoveRequire={vi.fn()}
        onReferenceOutput={vi.fn()}
      />,
    );
    fireEvent.click(screen.getByText(/Delete node/));
    expect(onDeleteNode).toHaveBeenCalledWith('alpha');
    expect(onClose).toHaveBeenCalled();
  });

  it('removes a requires edge and closes when the edge item is clicked', () => {
    const onRemoveRequire = vi.fn();
    const onClose = vi.fn();
    render(
      <WorkflowContextMenu
        menu={{ kind: WorkflowContextMenuKind.Edge, source: 'a', target: 'b', x: 10, y: 20 }}
        referenceOptions={[]}
        onClose={onClose}
        onDeleteNode={vi.fn()}
        onDeleteField={vi.fn()}
        onRemoveRequire={onRemoveRequire}
        onReferenceOutput={vi.fn()}
      />,
    );
    fireEvent.click(screen.getByText(/Remove requirement/));
    expect(onRemoveRequire).toHaveBeenCalledWith('a', 'b');
    expect(onClose).toHaveBeenCalled();
  });
});
