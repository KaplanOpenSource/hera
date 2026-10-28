import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, waitFor, cleanup, act } from '@testing-library/react';
import { mockApiCall } from './helpers/mockApiCall';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

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

let getJson: ReturnType<typeof vi.fn>;

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
  ({ getJson } = mockApiCall());
  useWorkflowRunStore.setState({ runs: {} });
});

describe('WorkflowRunPoller', () => {
  it('polls the running run and writes done + chunks back to the store', async () => {
    const chunks = [{ name: '__between__', text: 'the output' }];
    getJson.mockResolvedValueOnce({ status: 'done', error: '', chunks });

    render(<WorkflowRunPoller />);
    await act(async () => {
      startRun('w', 'tok');
    });

    await waitFor(() => {
      expect(getJson).toHaveBeenCalledWith('/workflow_status/tok', {}, { label: 'poll workflow', silent: true });
      expect(runOf('w').status).toBe(WorkflowRunStatus.Done);
      expect(runOf('w').chunks).toEqual(chunks);
      expect(mockPushInfo).toHaveBeenCalledWith('Workflow "w" finished');
    });
  });

  it('keeps polling while running, then stops on done', async () => {
    // Fake timers so the test does not sit through the real POLL_MS wait.
    vi.useFakeTimers();
    getJson.mockResolvedValueOnce({ status: 'running', error: '', chunks: [] });
    getJson.mockResolvedValueOnce({ status: 'done', error: '', chunks: [{ name: '__between__', text: 'done now' }] });

    render(<WorkflowRunPoller />);
    await act(async () => {
      startRun('w', 'tok');
    });

    // The first poll said "running", so the run is still going and a retry is pending.
    expect(runOf('w').status).toBe(WorkflowRunStatus.Running);
    expect(getJson).toHaveBeenCalledTimes(1);

    await act(async () => {
      await vi.advanceTimersByTimeAsync(POLL_MS);
    });

    expect(runOf('w').status).toBe(WorkflowRunStatus.Done);
    expect(getJson).toHaveBeenCalledTimes(2);
  });

  it('writes partial chunks to the store while still running', async () => {
    const partial = [{ name: '__preamble__', text: 'partial so far' }];
    const final = [{ name: '__between__', text: 'all done' }];
    getJson.mockResolvedValueOnce({ status: 'running', error: '', chunks: partial });
    getJson.mockResolvedValueOnce({ status: 'done', error: '', chunks: final });

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
    getJson.mockResolvedValueOnce({ status: 'error', error: 'it broke', chunks: [] });

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
    getJson.mockResolvedValueOnce({ status: 'not_found', error: '', chunks: null });

    render(<WorkflowRunPoller />);
    await act(async () => {
      startRun('w', 'gone');
    });

    await waitFor(() => {
      expect(runOf('w').status).toBe(WorkflowRunStatus.Error);
      expect(mockPushError).toHaveBeenCalledWith(expect.stringContaining('lost'));
    });
    expect(getJson).toHaveBeenCalledTimes(1);
  });

  it('marks the run as error when the poll call gets no answer from the server', async () => {
    getJson.mockResolvedValueOnce(null);

    render(<WorkflowRunPoller />);
    await act(async () => {
      startRun('w', 'tok');
    });

    await waitFor(() => {
      expect(runOf('w').status).toBe(WorkflowRunStatus.Error);
      expect(mockPushError).toHaveBeenCalledWith(expect.stringContaining('could not be read'));
    });
  });
});
