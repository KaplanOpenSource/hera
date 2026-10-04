import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, fireEvent, render, screen, within } from '@testing-library/react';
import { ReactFlow, ReactFlowProvider } from '@xyflow/react';
import { WorkflowNodeListPanel } from '../src/components/workflow/WorkflowNodeListPanel';

afterEach(() => cleanup());

const renderPanel = (nodeNames: string[], onPickNode = vi.fn(), handlers: {
  pinnedNodes?: string[],
  onTogglePin?: (name: string) => void,
  onDeleteNode?: (name: string) => void,
} = {}) => {
  render(
    <ReactFlowProvider>
      <ReactFlow nodes={[]} edges={[]}>
        <WorkflowNodeListPanel
          nodeNames={nodeNames}
          pinnedNodes={handlers.pinnedNodes ?? []}
          onPickNode={onPickNode}
          onTogglePin={handlers.onTogglePin ?? vi.fn()}
          onDeleteNode={handlers.onDeleteNode ?? vi.fn()}
        />
      </ReactFlow>
    </ReactFlowProvider>,
  );
  return onPickNode;
};

describe('WorkflowNodeListPanel', () => {
  it('starts collapsed and opens on the icon', () => {
    renderPanel(['a', 'b']);
    expect(screen.queryByText('a')).toBeNull();
    fireEvent.click(screen.getByLabelText('node list'));
    expect(screen.getByText('a')).toBeDefined();
    expect(screen.getByText('b')).toBeDefined();
  });

  it('reports the node clicked in the list', () => {
    const onPickNode = renderPanel(['a', 'b']);
    fireEvent.click(screen.getByLabelText('node list'));
    fireEvent.click(screen.getByText('b'));
    expect(onPickNode).toHaveBeenCalledWith('b');
  });

  it('puts the cursor in the search box when it opens', () => {
    renderPanel(['alpha']);
    fireEvent.click(screen.getByLabelText('node list'));
    expect(document.activeElement).toBe(screen.getByLabelText('search nodes'));
  });

  it('filters the list by the search text', () => {
    renderPanel(['alpha', 'beta']);
    fireEvent.click(screen.getByLabelText('node list'));
    fireEvent.change(screen.getByLabelText('search nodes'), { target: { value: 'BET' } });
    expect(screen.queryByText('alpha')).toBeNull();
    expect(screen.getByText('beta')).toBeDefined();
  });

  it('says so when nothing matches the search', () => {
    renderPanel(['alpha']);
    fireEvent.click(screen.getByLabelText('node list'));
    fireEvent.change(screen.getByLabelText('search nodes'), { target: { value: 'zz' } });
    expect(screen.getByText('No match.')).toBeDefined();
  });

  it('pins the node of the row from the list', () => {
    const onTogglePin = vi.fn();
    renderPanel(['alpha', 'beta'], vi.fn(), { onTogglePin });
    fireEvent.click(screen.getByLabelText('node list'));
    fireEvent.click(screen.getAllByLabelText('keep open')[1]);
    expect(onTogglePin).toHaveBeenCalledWith('beta');
  });

  it('deletes the node of the row from the list', () => {
    const onDeleteNode = vi.fn();
    const onPickNode = vi.fn();
    renderPanel(['alpha', 'beta'], onPickNode, { onDeleteNode });
    fireEvent.click(screen.getByLabelText('node list'));
    const rows = screen.getAllByRole('listitem');
    fireEvent.click(within(rows[0]).getByTestId('CloseIcon'));
    expect(onDeleteNode).toHaveBeenCalledWith('alpha');
    expect(onPickNode).not.toHaveBeenCalled();
  });

  it('says so when the workflow has no nodes', () => {
    renderPanel([]);
    fireEvent.click(screen.getByLabelText('node list'));
    expect(screen.getByText('No nodes.')).toBeDefined();
  });
});
