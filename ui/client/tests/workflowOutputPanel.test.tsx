import { describe, it, expect, afterEach, beforeEach } from 'vitest';
import { render, screen, cleanup, act } from '@testing-library/react';
import { WorkflowOutputPanel } from '../src/components/workflow/log/WorkflowOutputPanel';
import { useWorkflowRunStore, WorkflowRunStatus } from '../src/stores/useWorkflowRunStore';

const setRun = (status: WorkflowRunStatus, error: string, chunks: { name: string, text: string }[]) => {
  useWorkflowRunStore.setState({ runs: { w: { token: 't', status, error, chunks } } });
};

beforeEach(() => {
  useWorkflowRunStore.setState({ runs: {} });
});

afterEach(() => {
  cleanup();
});

describe('WorkflowOutputPanel', () => {
  it('says there is no run yet when the workflow never ran', () => {
    render(<WorkflowOutputPanel workflowName="w" />);

    expect(screen.getByText('No run yet.')).toBeTruthy();
  });

  it('shows the partial log and a running hint while running', () => {
    setRun(WorkflowRunStatus.Running, '', [{ name: '__preamble__', text: 'partial line while running' }]);
    render(<WorkflowOutputPanel workflowName="w" />);

    expect(screen.getByText('Running…')).toBeTruthy();
    expect(screen.getByText('partial line while running')).toBeTruthy();
  });

  it('shows the final log when done', () => {
    setRun(WorkflowRunStatus.Done, '', [{ name: '__between__', text: 'final output line' }]);
    render(<WorkflowOutputPanel workflowName="w" />);

    expect(screen.queryByText('Running…')).toBeNull();
    expect(screen.getByText('final output line')).toBeTruthy();
  });

  it('shows the error when the run failed', () => {
    setRun(WorkflowRunStatus.Error, 'it broke', []);
    render(<WorkflowOutputPanel workflowName="w" />);

    expect(screen.getByText('it broke')).toBeTruthy();
  });

  it('keeps the log visible alongside the error on failure', () => {
    setRun(WorkflowRunStatus.Error, 'it broke', [{ name: '__between__', text: 'log that led to the failure' }]);
    render(<WorkflowOutputPanel workflowName="w" />);

    expect(screen.getByText('it broke')).toBeTruthy();
    expect(screen.getByText('log that led to the failure')).toBeTruthy();
  });

  it('updates live as the store gains chunks', () => {
    setRun(WorkflowRunStatus.Running, '', [{ name: '__between__', text: 'first line' }]);
    render(<WorkflowOutputPanel workflowName="w" />);
    expect(screen.getByText('first line')).toBeTruthy();

    act(() => {
      useWorkflowRunStore.getState().setRunChunks('w', [
        { name: '__between__', text: 'first line' },
        { name: '__between__', text: 'second line' },
      ]);
    });

    expect(screen.getByText('second line')).toBeTruthy();
  });
});
