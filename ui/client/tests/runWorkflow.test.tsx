import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockStartWorkflow = vi.fn();
const mockPollWorkflow = vi.fn();
vi.mock('../src/io/runWorkflow', () => ({
  startWorkflow: (...args: any[]) => mockStartWorkflow(...args),
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

const { RunWorkflowButton } = await import('../src/components/workflow/RunWorkflowButton');
const { useViewSettingsStore } = await import('../src/stores/useViewSettingsStore');

const setSaving = (value: boolean) => {
  useViewSettingsStore.getState().setViewSettings({ alwaysSaveBeforeRun: value });
};
const isSaving = () => {
  return useViewSettingsStore.getState().viewSettings.alwaysSaveBeforeRun;
};

afterEach(() => {
  cleanup();
});

beforeEach(() => {
  vi.clearAllMocks();
  setSaving(false);
  // A run that starts and finishes on the first poll, unless a test overrides it.
  mockStartWorkflow.mockResolvedValue({ token: 't1' });
  mockPollWorkflow.mockResolvedValue({ status: 'done', output: '', error: '' });
});

describe('RunWorkflowButton', () => {
  it('starts the run, polls, and shows the output on done', async () => {
    mockStartWorkflow.mockResolvedValueOnce({ token: 'abc123' });
    mockPollWorkflow.mockResolvedValueOnce({ status: 'done', output: 'hello from hera\n', error: '' });

    render(<RunWorkflowButton projectName="TestProject" workflowName="hello_1" />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => {
      expect(mockStartWorkflow).toHaveBeenCalledWith({ projectName: 'TestProject', workflowName: 'hello_1' });
      expect(mockPollWorkflow).toHaveBeenCalledWith('abc123');
      expect(mockPushInfo).toHaveBeenCalledWith('Workflow "hello_1" finished');
      expect(screen.getByText(/hello from hera/i)).toBeTruthy();
    });
  });

  it('keeps polling while the run is still running, then stops on done', async () => {
    mockStartWorkflow.mockResolvedValueOnce({ token: 'tok' });
    mockPollWorkflow.mockResolvedValueOnce({ status: 'running', output: '', error: '' });
    mockPollWorkflow.mockResolvedValueOnce({ status: 'done', output: 'finished output', error: '' });

    render(<RunWorkflowButton projectName="P" workflowName="w" />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => {
      expect(screen.getByText(/finished output/i)).toBeTruthy();
    });
    // Exactly the running poll + the done poll; it does not keep polling after done.
    expect(mockPollWorkflow).toHaveBeenCalledTimes(2);
  });

  it('shows an error when a poll reports the run failed', async () => {
    mockStartWorkflow.mockResolvedValueOnce({ token: 'tok' });
    mockPollWorkflow.mockResolvedValueOnce({ status: 'error', output: '', error: 'boom in the run' });

    render(<RunWorkflowButton projectName="P" workflowName="w" />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledWith('run workflow: boom in the run');
      expect(screen.getByText(/boom in the run/i)).toBeTruthy();
    });
    expect(mockPushInfo).not.toHaveBeenCalled();
  });

  it('shows a busy message and does not poll when the server is busy', async () => {
    mockStartWorkflow.mockResolvedValueOnce({ status: 'busy' });

    render(<RunWorkflowButton projectName="P" workflowName="w" />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledWith(expect.stringContaining('busy'));
    });
    expect(mockPollWorkflow).not.toHaveBeenCalled();
    expect(mockPushInfo).not.toHaveBeenCalled();
  });

  it('stops polling on an unknown token', async () => {
    mockStartWorkflow.mockResolvedValueOnce({ token: 'gone' });
    mockPollWorkflow.mockResolvedValueOnce({ status: 'not_found', output: '', error: '' });

    render(<RunWorkflowButton projectName="P" workflowName="w" />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledWith(expect.stringContaining('lost'));
    });
    // One poll only: it stops instead of looping on a token the server does not know.
    expect(mockPollWorkflow).toHaveBeenCalledTimes(1);
    expect(mockPushInfo).not.toHaveBeenCalled();
  });

  it('shows an error snackbar when starting the run fails', async () => {
    mockStartWorkflow.mockReset();
    mockStartWorkflow.mockRejectedValueOnce(new Error('boom'));

    render(<RunWorkflowButton projectName="TestProject" workflowName="hello_1" />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledWith('run workflow: boom');
    });
    expect(mockPollWorkflow).not.toHaveBeenCalled();
    expect(mockPushInfo).not.toHaveBeenCalled();
  });

  it('is disabled with a reason when told to be', () => {
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

  it('does not save on a plain left click when there are no unsaved changes', async () => {
    const save = vi.fn().mockResolvedValue(undefined);

    render(<RunWorkflowButton projectName="P" workflowName="w" save={save} />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => expect(mockStartWorkflow).toHaveBeenCalled());
    expect(save).not.toHaveBeenCalled();
  });

  it('disables the run button when there are unsaved changes and saving is off', () => {
    const save = vi.fn().mockResolvedValue(undefined);

    render(<RunWorkflowButton projectName="P" workflowName="w" isChanged save={save} />);
    expect(screen.getByRole('button', { name: /save changes before running/i })).toHaveProperty('disabled', true);
  });

  it('enables the run button with unsaved changes once saving is on', () => {
    setSaving(true);
    const save = vi.fn().mockResolvedValue(undefined);

    render(<RunWorkflowButton projectName="P" workflowName="w" isChanged save={save} />);
    expect(screen.getByRole('button', { name: /run workflow/i })).toHaveProperty('disabled', false);
  });

  it('saves before running on left click when the flag is on', async () => {
    setSaving(true);
    const order: string[] = [];
    const save = vi.fn().mockImplementation(async () => { order.push('save'); });
    mockStartWorkflow.mockReset();
    mockStartWorkflow.mockImplementation(async () => { order.push('run'); return { token: 't' }; });

    render(<RunWorkflowButton projectName="P" workflowName="w" isChanged save={save} />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => expect(mockStartWorkflow).toHaveBeenCalled());
    expect(order).toEqual(['save', 'run']);
  });

  it('runs with save from the right-click menu (one time, flag stays off)', async () => {
    const save = vi.fn().mockResolvedValue(undefined);

    render(<RunWorkflowButton projectName="P" workflowName="w" isChanged save={save} />);
    fireEvent.contextMenu(screen.getByRole('button'));
    await act(async () => {
      fireEvent.click(screen.getByText('Run with save'));
    });

    await waitFor(() => expect(save).toHaveBeenCalledTimes(1));
    expect(mockStartWorkflow).toHaveBeenCalled();
    expect(isSaving()).toBe(false);
  });

  it('toggles the always-save flag on and off from the menu', async () => {
    const save = vi.fn().mockResolvedValue(undefined);

    render(<RunWorkflowButton projectName="P" workflowName="w" isChanged save={save} />);

    // Turn it on.
    fireEvent.contextMenu(screen.getByRole('button'));
    fireEvent.click(screen.getByText('Auto save before run'));
    expect(isSaving()).toBe(true);

    // The toggle now offers the opposite action.
    fireEvent.contextMenu(screen.getByRole('button'));
    expect(screen.queryByText('Auto save before run')).toBeNull();
    fireEvent.click(screen.getByText('Stop auto save before run'));
    expect(isSaving()).toBe(false);
  });

  it('keeps every run button in sync when the flag toggles', () => {
    const save = vi.fn().mockResolvedValue(undefined);

    render(
      <>
        <RunWorkflowButton projectName="P" workflowName="w" isChanged save={save} />
        <RunWorkflowButton projectName="P" workflowName="w" isChanged save={save} />
      </>,
    );
    const buttons = () => screen.getAllByRole('button') as HTMLButtonElement[];
    // Both start disabled: unsaved changes with saving off.
    expect(buttons()).toHaveLength(2);
    expect(buttons().every(b => b.disabled)).toBe(true);

    // Turn saving on from the first button's menu.
    fireEvent.contextMenu(buttons()[0]);
    fireEvent.click(screen.getByText('Auto save before run'));

    // Both buttons react, not just the one that toggled it.
    expect(buttons().every(b => !b.disabled)).toBe(true);
  });

  it('disables "Run with save" when there are no unsaved changes', () => {
    const save = vi.fn().mockResolvedValue(undefined);

    render(<RunWorkflowButton projectName="P" workflowName="w" save={save} />);
    fireEvent.contextMenu(screen.getByRole('button'));
    expect(screen.getByText('Run with save').closest('li')?.getAttribute('aria-disabled')).toBe('true');
  });
});
