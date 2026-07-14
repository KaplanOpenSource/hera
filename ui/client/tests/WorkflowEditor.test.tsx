import { describe, it, expect, vi, beforeEach } from 'vitest';
import { act, cleanup, render, screen } from '@testing-library/react';

// The editor reads the node catalog from a zustand store and mounts a reader
// that fetches it; stub both so nothing hits the network.
vi.mock('../src/components/workflow/useNodeCatalog', () => ({
  useNodeCatalog: (selector: (state: { catalog: unknown[] }) => unknown) => selector({ catalog: [] }),
  NodeCatalogReader: () => null,
}));

// Replace the ReactFlow-based graph with a stub that just captures the handler
// props, so we can drive WorkflowEditor's logic without rendering ReactFlow.
let graphProps: any = null;
vi.mock('../src/components/workflow/WorkflowGraph', () => ({
  WorkflowGraph: (props: any) => {
    graphProps = props;
    return null;
  },
}));

import { WorkflowEditor } from '../src/components/workflow/WorkflowEditor';

const setup = (workflow: any) => {
  const setWorkflow = vi.fn();
  render(<WorkflowEditor workflow={workflow} setWorkflow={setWorkflow} />);
  return { setWorkflow };
};

beforeEach(() => {
  graphProps = null;
  cleanup();
});

describe('WorkflowEditor empty state', () => {
  it('shows a message and no graph when there is no workflow block', () => {
    setup({});
    expect(screen.getByText('No workflow found in this document.')).toBeDefined();
    expect(graphProps).toBeNull();
  });
});

describe('WorkflowEditor addNode', () => {
  it('adds node1 with empty defaults to an empty block', () => {
    const { setWorkflow } = setup({ nodes: {}, nodeList: [] });
    act(() => graphProps.onAddNode());
    expect(setWorkflow).toHaveBeenCalledWith({
      nodeList: ['node1'],
      nodes: { node1: { type: '', Execution: { input_parameters: {} } } },
    });
  });

  it('picks the next free name when node1 already exists', () => {
    const { setWorkflow } = setup({ nodeList: ['node1'], nodes: { node1: {} } });
    act(() => graphProps.onAddNode());
    const written = setWorkflow.mock.calls[0][0];
    expect(written.nodeList).toEqual(['node1', 'node2']);
    expect(written.nodes.node2).toEqual({ type: '', Execution: { input_parameters: {} } });
  });

  it('preserves an extra { workflow } nesting level', () => {
    const { setWorkflow } = setup({ workflow: { nodeList: ['a'], nodes: { a: {} } } });
    act(() => graphProps.onAddNode());
    const written = setWorkflow.mock.calls[0][0];
    expect(written.workflow.nodeList).toEqual(['a', 'node2']);
    expect(written.workflow.nodes.node2).toBeDefined();
  });
});

describe('WorkflowEditor renameNode', () => {
  it('renames the key, the node list, and references in requires', () => {
    const { setWorkflow } = setup({
      nodeList: ['a', 'b'],
      nodes: { a: { type: 't' }, b: { type: 't2', requires: 'a' } },
    });
    act(() => graphProps.onRenameNode('a', 'x'));
    expect(setWorkflow).toHaveBeenCalledWith({
      nodeList: ['x', 'b'],
      nodes: { x: { type: 't' }, b: { type: 't2', requires: 'x' } },
    });
  });

  it('rewrites a requires array', () => {
    const { setWorkflow } = setup({
      nodeList: ['a', 'b'],
      nodes: { a: {}, b: { requires: ['a', 'c'] } },
    });
    act(() => graphProps.onRenameNode('a', 'x'));
    expect(setWorkflow.mock.calls[0][0].nodes.b.requires).toEqual(['x', 'c']);
  });

  it('does nothing for an empty, unchanged, or duplicate name', () => {
    const { setWorkflow } = setup({ nodeList: ['a', 'b'], nodes: { a: {}, b: {} } });
    act(() => graphProps.onRenameNode('a', '   '));
    act(() => graphProps.onRenameNode('a', 'a'));
    act(() => graphProps.onRenameNode('a', 'b'));
    expect(setWorkflow).not.toHaveBeenCalled();
  });
});

describe('WorkflowEditor deleteNode', () => {
  it('removes the node from the list and the map', () => {
    const { setWorkflow } = setup({ nodeList: ['a', 'b'], nodes: { a: {}, b: {} } });
    act(() => graphProps.onDeleteNode('a'));
    expect(setWorkflow).toHaveBeenCalledWith({ nodeList: ['b'], nodes: { b: {} } });
  });
});

describe('WorkflowEditor requires edges', () => {
  it('adds a requires edge to the target node', () => {
    const { setWorkflow } = setup({ nodeList: ['a', 'b'], nodes: { a: {}, b: {} } });
    act(() => graphProps.onAddRequire('a', 'b'));
    expect(setWorkflow.mock.calls[0][0].nodes.b.requires).toEqual(['a']);
  });

  it('does not add a duplicate requires edge', () => {
    const { setWorkflow } = setup({ nodeList: ['a', 'b'], nodes: { a: {}, b: { requires: ['a'] } } });
    act(() => graphProps.onAddRequire('a', 'b'));
    expect(setWorkflow).not.toHaveBeenCalled();
  });

  it('drops the requires field when the last edge is removed', () => {
    const { setWorkflow } = setup({ nodeList: ['a', 'b'], nodes: { a: {}, b: { requires: ['a'] } } });
    act(() => graphProps.onRemoveRequire('a', 'b'));
    expect(setWorkflow.mock.calls[0][0].nodes.b).toEqual({});
  });

  it('keeps the remaining edges when one is removed', () => {
    const { setWorkflow } = setup({ nodeList: ['a', 'b'], nodes: { a: {}, b: { requires: ['a', 'c'] } } });
    act(() => graphProps.onRemoveRequire('a', 'b'));
    expect(setWorkflow.mock.calls[0][0].nodes.b.requires).toEqual(['c']);
  });
});

describe('WorkflowEditor setNode', () => {
  it('replaces a node and keeps the rest of the block', () => {
    const { setWorkflow } = setup({ nodeList: ['a'], nodes: { a: { type: 'old' } } });
    act(() => graphProps.onSetNode('a', { type: 'new' }));
    expect(setWorkflow).toHaveBeenCalledWith({ nodeList: ['a'], nodes: { a: { type: 'new' } } });
  });
});
