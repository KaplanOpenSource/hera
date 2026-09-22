import { describe, it, expect, beforeEach } from 'vitest';
import { useWorkflowFocusStore } from '../src/stores/useWorkflowFocusStore';

const focusNode = (workflowName: string, nodeName: string) => {
  useWorkflowFocusStore.getState().focusNode(workflowName, nodeName);
};

const focus = () => {
  return useWorkflowFocusStore.getState().focus;
};

beforeEach(() => {
  useWorkflowFocusStore.setState({ focus: null });
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
