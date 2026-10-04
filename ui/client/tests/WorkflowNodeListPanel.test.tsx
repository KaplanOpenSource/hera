import { describe, it, expect, vi, afterEach } from 'vitest';
import { cleanup, fireEvent, render, screen } from '@testing-library/react';
import { ReactFlow, ReactFlowProvider } from '@xyflow/react';
import { WorkflowNodeListPanel } from '../src/components/workflow/WorkflowNodeListPanel';

afterEach(() => cleanup());

const renderPanel = (nodeNames: string[], onPickNode = vi.fn()) => {
  render(
    <ReactFlowProvider>
      <ReactFlow nodes={[]} edges={[]}>
        <WorkflowNodeListPanel nodeNames={nodeNames} onPickNode={onPickNode} />
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

  it('says so when the workflow has no nodes', () => {
    renderPanel([]);
    fireEvent.click(screen.getByLabelText('node list'));
    expect(screen.getByText('No nodes.')).toBeDefined();
  });
});
