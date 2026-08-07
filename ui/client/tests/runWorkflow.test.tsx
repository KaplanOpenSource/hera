import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockRunWorkflow = vi.fn();
vi.mock('../src/io/runWorkflow', () => ({
  runWorkflow: (...args: any[]) => mockRunWorkflow(...args),
}));

const mockPushInfo = vi.fn();
const mockPushError = vi.fn();
vi.mock('../src/io/snackbar', () => ({
  pushInfo: (...args: any[]) => mockPushInfo(...args),
  pushError: (...args: any[]) => mockPushError(...args),
  pushRunning: () => 'key',
  dismiss: vi.fn(),
}));

const { RunWorkflowButton } = await import('../src/components/workflow/RunWorkflowButton');

afterEach(() => {
  cleanup();
});

beforeEach(() => {
  vi.clearAllMocks();
});

describe('RunWorkflowButton', () => {
  it('runs the workflow via the endpoint and shows its output', async () => {
    mockRunWorkflow.mockResolvedValueOnce({ dispatch_id: 'abc123', output: 'hello from hera\n' });

    render(<RunWorkflowButton projectName="TestProject" workflowName="hello_1" />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => {
      expect(mockRunWorkflow).toHaveBeenCalledWith({ projectName: 'TestProject', workflowName: 'hello_1' });
      expect(mockPushInfo).toHaveBeenCalledWith('Workflow "hello_1" finished');
      expect(screen.getByText(/hello from hera/i)).toBeTruthy();
    });
  });

  it('shows an error snackbar when the run fails', async () => {
    mockRunWorkflow.mockRejectedValueOnce(new Error('boom'));

    render(<RunWorkflowButton projectName="TestProject" workflowName="hello_1" />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledWith('run workflow: boom');
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
