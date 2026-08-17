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
    mockRunWorkflow.mockResolvedValueOnce({ dispatch_id: 'x', output: '' });
    const save = vi.fn().mockResolvedValue(undefined);

    render(<RunWorkflowButton projectName="P" workflowName="w" save={save} />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => expect(mockRunWorkflow).toHaveBeenCalled());
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
    mockRunWorkflow.mockImplementation(async () => { order.push('run'); return { dispatch_id: 'x', output: '' }; });

    render(<RunWorkflowButton projectName="P" workflowName="w" isChanged save={save} />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /run workflow/i }));
    });

    await waitFor(() => expect(mockRunWorkflow).toHaveBeenCalled());
    expect(order).toEqual(['save', 'run']);
  });

  it('runs with save from the right-click menu (one time, flag stays off)', async () => {
    mockRunWorkflow.mockResolvedValue({ dispatch_id: 'x', output: '' });
    const save = vi.fn().mockResolvedValue(undefined);

    render(<RunWorkflowButton projectName="P" workflowName="w" isChanged save={save} />);
    fireEvent.contextMenu(screen.getByRole('button'));
    await act(async () => {
      fireEvent.click(screen.getByText('Run with save'));
    });

    await waitFor(() => expect(save).toHaveBeenCalledTimes(1));
    expect(mockRunWorkflow).toHaveBeenCalled();
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
