import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetchPython = vi.fn();
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: (...args: any[]) => mockFetchPython(...args),
}));

const mockPushInfo = vi.fn();
vi.mock('../src/io/snackbar', () => ({
  pushInfo: (...args: any[]) => mockPushInfo(...args),
}));

const { RunWorkflowButton } = await import('../src/components/workflow/RunWorkflowButton');

afterEach(() => {
  cleanup();
});

beforeEach(() => {
  vi.clearAllMocks();
});

describe('RunWorkflowButton', () => {
  it('builds and executes the workflow from the DB with the local scheduler', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: { dispatch_id: 'abc123' }, problem: undefined });

    render(<RunWorkflowButton projectName="TestProject" workflowName="hello_1" />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
      const code = mockFetchPython.mock.calls[0][0].code;
      expect(code).toContain('SIMULATIONS_WORKFLOWS');
      expect(code).toContain("projectName='TestProject'");
      expect(code).toContain("executeWorkflowFromDB('hello_1', scheduler='local')");
      expect(mockPushInfo).toHaveBeenCalledWith('Workflow "hello_1" finished');
    });
  });

  it('does not report success when the run fails', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: undefined, problem: 'exec failed' });

    render(<RunWorkflowButton projectName="TestProject" workflowName="hello_1" />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });
    expect(mockPushInfo).not.toHaveBeenCalled();
  });

  it('is disabled with a reason when there are unsaved changes', () => {
    render(
      <RunWorkflowButton
        projectName="TestProject"
        workflowName="hello_1"
        disabled
        disabledReason="Save changes before running"
      />,
    );
    expect(screen.getByRole('button', { name: /save changes before running/i })).toHaveProperty('disabled', true);
  });
});
