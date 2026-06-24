import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, fireEvent, render, screen } from '@testing-library/react';
import { WorkflowContextMenu, WorkflowContextMenuKind } from '../src/components/workflow/WorkflowContextMenu';

afterEach(() => cleanup());

describe('WorkflowContextMenu', () => {
  it('renders nothing when there is no target', () => {
    render(<WorkflowContextMenu menu={null} onClose={vi.fn()} onDeleteNode={vi.fn()} onRemoveRequire={vi.fn()} />);
    expect(screen.queryByText(/Delete node/)).toBeNull();
    expect(screen.queryByText(/Remove requirement/)).toBeNull();
  });

  it('deletes a node and closes when the node item is clicked', () => {
    const onDeleteNode = vi.fn();
    const onClose = vi.fn();
    render(
      <WorkflowContextMenu
        menu={{ kind: WorkflowContextMenuKind.Node, name: 'alpha', x: 10, y: 20 }}
        onClose={onClose}
        onDeleteNode={onDeleteNode}
        onRemoveRequire={vi.fn()}
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
        onClose={onClose}
        onDeleteNode={vi.fn()}
        onRemoveRequire={onRemoveRequire}
      />,
    );
    fireEvent.click(screen.getByText(/Remove requirement/));
    expect(onRemoveRequire).toHaveBeenCalledWith('a', 'b');
    expect(onClose).toHaveBeenCalled();
  });
});
