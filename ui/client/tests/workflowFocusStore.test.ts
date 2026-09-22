import { describe, it, expect, beforeEach } from 'vitest';
import { useWorkflowFocusStore } from '../src/stores/useWorkflowFocusStore';

const focusNode = (workflowName: string, nodeName: string) => {
  useWorkflowFocusStore.getState().focusNode(workflowName, nodeName);
};

const focus = () => {
  return useWorkflowFocusStore.getState().focus;
};

beforeEach(() => {
  useWorkflowFocusStore.setState({ focus: null, hover: null });
});

describe('useWorkflowFocusStore', () => {
  it('starts with nothing focused', () => {
    expect(focus()).toBeNull();
  });

  it('records the workflow and node that was picked', () => {
    focusNode('Workflow2', 'ListFiles');

    expect(focus()?.workflowName).toBe('Workflow2');
    expect(focus()?.nodeName).toBe('ListFiles');
  });

  it('makes a repeat pick of the same node a fresh request', () => {
    focusNode('Workflow2', 'ListFiles');
    const first = focus()!.seq;
    focusNode('Workflow2', 'ListFiles');

    expect(focus()!.seq).toBeGreaterThan(first);
  });

  it('replaces the request when another workflow picks a node', () => {
    focusNode('Workflow2', 'ListFiles');
    focusNode('Other', 'Build');

    expect(focus()?.workflowName).toBe('Other');
    expect(focus()?.nodeName).toBe('Build');
  });
});

describe('useWorkflowFocusStore hover', () => {
  const hoverNode = (workflowName: string, nodeName: string | null) => {
    useWorkflowFocusStore.getState().hoverNode(workflowName, nodeName);
  };
  const hover = () => {
    return useWorkflowFocusStore.getState().hover;
  };

  it('starts with nothing hovered', () => {
    expect(hover()).toBeNull();
  });

  it('records the hovered node', () => {
    hoverNode('Workflow2', 'ListFiles');

    expect(hover()).toEqual({ workflowName: 'Workflow2', nodeName: 'ListFiles' });
  });

  it('clears when the pointer leaves the node', () => {
    hoverNode('Workflow2', 'ListFiles');
    hoverNode('Workflow2', null);

    expect(hover()).toBeNull();
  });
});
