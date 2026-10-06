import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, fireEvent, render, screen } from '@testing-library/react';
import { WorkflowContextMenu, WorkflowContextMenuKind } from '../src/components/workflow/WorkflowContextMenu';
import { Reference } from '../src/components/workflow/references/Reference';
import { INPUT, OUTPUT } from '../src/components/workflow/references/knownKinds';

afterEach(() => cleanup());

const referenceOptions = [
  { node: 'src', references: [new Reference('src', OUTPUT, 'out1'), new Reference('src', OUTPUT, 'out2')] },
  { node: 'other', references: [new Reference('other', OUTPUT, 'bee')] },
  { node: 'third', references: [new Reference('third', INPUT, 'cmd')] },
];

describe('WorkflowContextMenu', () => {
  it('renders nothing when there is no target', () => {
    render(<WorkflowContextMenu menu={null} referenceOptions={[]} onClose={vi.fn()} onDeleteNode={vi.fn()} onDeleteField={vi.fn()} onRemoveRequire={vi.fn()} onReferenceOutput={vi.fn()} />);
    expect(screen.queryByText(/Delete node/)).toBeNull();
    expect(screen.queryByText(/Remove requirement/)).toBeNull();
  });

  it('references another node output via two fly-out submenus', async () => {
    const onReferenceOutput = vi.fn();
    const onClose = vi.fn();
    render(
      <WorkflowContextMenu
        menu={{ kind: WorkflowContextMenuKind.Field, node: 'alpha', param: 'p', x: 10, y: 20, caret: 3 }}
        referenceOptions={referenceOptions}
        onClose={onClose}
        onDeleteNode={vi.fn()}
        onDeleteField={vi.fn()}
        onRemoveRequire={vi.fn()}
        onReferenceOutput={onReferenceOutput}
      />,
    );
    // Neither submenu autocomplete shows until the item is opened.
    expect(screen.queryByRole('combobox', { name: 'Node' })).toBeNull();

    // Open the "Reference another node" fly-out → the node autocomplete appears.
    fireEvent.click(screen.getByText('Reference another node'));
    expect(screen.queryByRole('combobox', { name: 'Parameter' })).toBeNull();

    // Pick the source node → the output submenu flies out.
    fireEvent.change(screen.getByRole('combobox', { name: 'Node' }), { target: { value: 'src' } });
    fireEvent.click(await screen.findByRole('option', { name: 'src' }));

    // Pick one of that node's outputs → inserts the reference and closes.
    fireEvent.change(await screen.findByRole('combobox', { name: 'Parameter' }), { target: { value: 'out2' } });
    fireEvent.click(await screen.findByRole('option', { name: 'out2' }));
    expect(onReferenceOutput).toHaveBeenCalledWith('alpha', 'p', new Reference('src', OUTPUT, 'out2'), 3);
    expect(onClose).toHaveBeenCalled();
  });

  it('references another node input parameter', async () => {
    const onReferenceOutput = vi.fn();
    render(
      <WorkflowContextMenu
        menu={{ kind: WorkflowContextMenuKind.Field, node: 'alpha', param: 'p', x: 10, y: 20, caret: 3 }}
        referenceOptions={referenceOptions}
        onClose={vi.fn()}
        onDeleteNode={vi.fn()}
        onDeleteField={vi.fn()}
        onRemoveRequire={vi.fn()}
        onReferenceOutput={onReferenceOutput}
      />,
    );
    fireEvent.click(screen.getByText('Reference another node'));
    fireEvent.change(screen.getByRole('combobox', { name: 'Node' }), { target: { value: 'third' } });
    fireEvent.click(await screen.findByRole('option', { name: 'third' }));
    fireEvent.change(await screen.findByRole('combobox', { name: 'Parameter' }), { target: { value: 'cmd' } });
    fireEvent.click(await screen.findByRole('option', { name: 'cmd' }));
    expect(onReferenceOutput).toHaveBeenCalledWith('alpha', 'p', new Reference('third', INPUT, 'cmd'), 3);
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
