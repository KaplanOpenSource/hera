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

const { WorkflowRunPoller, POLL_MS } = await import('../src/components/workflow/WorkflowRunPoller');
const { useWorkflowRunStore, WorkflowRunStatus } = await import('../src/stores/useWorkflowRunStore');

const startRun = (workflowName: string, token: string) => {
  useWorkflowRunStore.getState().startRun(workflowName, token);
};
const runOf = (workflowName: string) => {
  return useWorkflowRunStore.getState().runs[workflowName];
};

afterEach(() => {
  cleanup();
  vi.useRealTimers();
});

beforeEach(() => {
  vi.clearAllMocks();
  useWorkflowRunStore.setState({ runs: {} });
});

describe('WorkflowRunPoller', () => {
  it('polls the running run and writes done + chunks back to the store', async () => {
    const chunks = [{ name: '__between__', text: 'the output' }];
    mockPollWorkflow.mockResolvedValueOnce({ status: 'done', error: '', chunks });

    render(<WorkflowRunPoller />);
    await act(async () => {
      startRun('w', 'tok');
    });

    await waitFor(() => {
      expect(mockPollWorkflow).toHaveBeenCalledWith('tok');
      expect(runOf('w').status).toBe(WorkflowRunStatus.Done);
      expect(runOf('w').chunks).toEqual(chunks);
      expect(mockPushInfo).toHaveBeenCalledWith('Workflow "w" finished');
    });
  });

  it('keeps polling while running, then stops on done', async () => {
    // Fake timers so the test does not sit through the real POLL_MS wait.
    vi.useFakeTimers();
    mockPollWorkflow.mockResolvedValueOnce({ status: 'running', error: '', chunks: [] });
    mockPollWorkflow.mockResolvedValueOnce({ status: 'done', error: '', chunks: [{ name: '__between__', text: 'done now' }] });

    render(<WorkflowRunPoller />);
    await act(async () => {
      startRun('w', 'tok');
    });

    // The first poll said "running", so the run is still going and a retry is pending.
    expect(runOf('w').status).toBe(WorkflowRunStatus.Running);
    expect(mockPollWorkflow).toHaveBeenCalledTimes(1);

    await act(async () => {
      await vi.advanceTimersByTimeAsync(POLL_MS);
    });

    expect(runOf('w').status).toBe(WorkflowRunStatus.Done);
    expect(mockPollWorkflow).toHaveBeenCalledTimes(2);
  });

  it('writes partial chunks to the store while still running', async () => {
    const partial = [{ name: '__preamble__', text: 'partial so far' }];
    const final = [{ name: '__between__', text: 'all done' }];
    mockPollWorkflow.mockResolvedValueOnce({ status: 'running', error: '', chunks: partial });
    mockPollWorkflow.mockResolvedValueOnce({ status: 'done', error: '', chunks: final });

    // Fake timers so the test does not sit through the real POLL_MS wait.
    vi.useFakeTimers();

    render(<WorkflowRunPoller />);
    await act(async () => {
      startRun('w', 'tok');
    });

    // The running poll updates chunks but keeps the run in the running state.
    expect(runOf('w').chunks).toEqual(partial);
    expect(runOf('w').status).toBe(WorkflowRunStatus.Running);

    await act(async () => {
      await vi.advanceTimersByTimeAsync(POLL_MS);
    });

    // Then it finishes with the final chunks.
    expect(runOf('w').status).toBe(WorkflowRunStatus.Done);
    expect(runOf('w').chunks).toEqual(final);
  });

  it('writes error to the store when a poll reports failure', async () => {
    mockPollWorkflow.mockResolvedValueOnce({ status: 'error', error: 'it broke', chunks: [] });

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
    mockPollWorkflow.mockResolvedValueOnce({ status: 'not_found', error: '', chunks: null });

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
