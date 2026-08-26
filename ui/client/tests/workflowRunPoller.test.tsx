import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockPollWorkflow = vi.fn();
vi.mock('../src/io/runWorkflow', () => ({
  pollWorkflow: (...args: any[]) => mockPollWorkflow(...args),
}));

const mockPushInfo = vi.fn();
const mockPushError = vi.fn();
vi.mock('../src/io/snackbar', () => ({
  pushInfo: (...args: any[]) => mockPushInfo(...args),
  pushError: (...args: any[]) => mockPushError(...args),
  pushRunning: () => 'key',
  dismiss: vi.fn(),
}));

const { WorkflowRunPoller } = await import('../src/components/workflow/WorkflowRunPoller');
const { useWorkflowRunStore, WorkflowRunStatus } = await import('../src/stores/useWorkflowRunStore');

const startRun = (workflowName: string, token: string) => {
  useWorkflowRunStore.getState().startRun(workflowName, token);
};
const runOf = (workflowName: string) => {
  return useWorkflowRunStore.getState().runs[workflowName];
};

afterEach(() => {
  cleanup();
});

beforeEach(() => {
  vi.clearAllMocks();
  useWorkflowRunStore.setState({ runs: {} });
});

describe('WorkflowRunPoller', () => {
  it('polls the running run and writes done + output back to the store', async () => {
    mockPollWorkflow.mockResolvedValueOnce({ status: 'done', output: 'the output', error: '' });

    render(<WorkflowRunPoller />);
    await act(async () => {
      startRun('w', 'tok');
    });

    await waitFor(() => {
      expect(mockPollWorkflow).toHaveBeenCalledWith('tok');
      expect(runOf('w').status).toBe(WorkflowRunStatus.Done);
      expect(runOf('w').output).toBe('the output');
      expect(mockPushInfo).toHaveBeenCalledWith('Workflow "w" finished');
    });
  });

  it('keeps polling while running, then stops on done', async () => {
    mockPollWorkflow.mockResolvedValueOnce({ status: 'running', output: '', error: '' });
    mockPollWorkflow.mockResolvedValueOnce({ status: 'done', output: 'done now', error: '' });

    render(<WorkflowRunPoller />);
    await act(async () => {
      startRun('w', 'tok');
    });

    await waitFor(() => expect(runOf('w').status).toBe(WorkflowRunStatus.Done));
    expect(mockPollWorkflow).toHaveBeenCalledTimes(2);
  });

  it('writes partial output to the store while still running', async () => {
    mockPollWorkflow.mockResolvedValueOnce({ status: 'running', output: 'partial so far', error: '' });
    mockPollWorkflow.mockResolvedValueOnce({ status: 'done', output: 'all done', error: '' });

    render(<WorkflowRunPoller />);
    await act(async () => {
      startRun('w', 'tok');
    });

    // The running poll updates output but keeps the run in the running state.
    await waitFor(() => expect(runOf('w').output).toBe('partial so far'));
    // Then it finishes with the final output.
    await waitFor(() => expect(runOf('w').status).toBe(WorkflowRunStatus.Done));
    expect(runOf('w').output).toBe('all done');
  });

  it('writes error to the store when a poll reports failure', async () => {
    mockPollWorkflow.mockResolvedValueOnce({ status: 'error', output: '', error: 'it broke' });

    render(<WorkflowRunPoller />);
    await act(async () => {
      startRun('w', 'tok');
    });

    await waitFor(() => {
      expect(runOf('w').status).toBe(WorkflowRunStatus.Error);
      expect(runOf('w').error).toBe('it broke');
      expect(mockPushError).toHaveBeenCalledWith('run workflow: it broke');
    });
    expect(mockPushInfo).not.toHaveBeenCalled();
  });

  it('stops on an unknown token and marks the run as error', async () => {
    mockPollWorkflow.mockResolvedValueOnce({ status: 'not_found', output: '', error: '' });

    render(<WorkflowRunPoller />);
    await act(async () => {
      startRun('w', 'gone');
    });

    await waitFor(() => {
      expect(runOf('w').status).toBe(WorkflowRunStatus.Error);
      expect(mockPushError).toHaveBeenCalledWith(expect.stringContaining('lost'));
    });
    expect(mockPollWorkflow).toHaveBeenCalledTimes(1);
  });
});
