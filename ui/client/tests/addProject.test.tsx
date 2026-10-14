import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

// The dialog also lists the registered repositories, which still uses /exec.
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: vi.fn().mockResolvedValue({ data: null }),
  fetchPythonClean: vi.fn().mockResolvedValue({ data: null }),
  WARMING_UP: 'WARMING_UP',
}));

const mockPostJson = vi.fn();
vi.mock('../src/io/entrypoint', () => ({
  postJson: (...args: any[]) => mockPostJson(...args),
  getJson: vi.fn(),
}));

const { AddProjectButton } = await import('../src/components/header/AddProjectButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');

afterEach(() => {
  cleanup();
});

beforeEach(() => {
  vi.clearAllMocks();
  mockPostJson.mockResolvedValue({
    project: { name: 'MyProject', documents: [] },
    addedNotebooks: 0,
  });
});

const openDialogAndSubmit = async (projectName: string, filesDirectory?: string) => {
  render(
    <MemoryRouter>
      <AddProjectButton />
    </MemoryRouter>
  );
  // The trigger is an icon-only button (no accessible name); it's the only button before the dialog opens.
  await act(async () => {
    fireEvent.click(screen.getAllByRole('button')[0]);
  });
  const dialog = await screen.findByRole('dialog');
  fireEvent.change(within(dialog).getByRole('textbox', { name: /project name/i }), { target: { value: projectName } });
  if (filesDirectory !== undefined) {
    fireEvent.change(within(dialog).getByRole('textbox', { name: /project files directory/i }), { target: { value: filesDirectory } });
  }
  await act(async () => {
    fireEvent.click(within(dialog).getByRole('button', { name: /add project/i }));
  });
};

describe('AddProjectButton', () => {
  it('creates the project in one call to project.create', async () => {
    await openDialogAndSubmit('MyProject');

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    const [path, body] = mockPostJson.mock.calls[0];
    expect(path).toBe('/project/create');
    expect(body).toEqual({
      projectName: 'MyProject',
      filesDirectory: '',
      loadRepositories: true,
    });
  });

  it('sends the files directory the user typed', async () => {
    await openDialogAndSubmit('MyProject', '/data/existing');

    await waitFor(() => {
      expect(mockPostJson.mock.calls[0][1].filesDirectory).toBe('/data/existing');
    });
  });

  it('shows the project the server answered with', async () => {
    const project = { name: 'MyProject', documents: [{ _id: { $oid: 'nb' } }] };
    mockPostJson.mockResolvedValueOnce({ project, addedNotebooks: 1 });

    await openDialogAndSubmit('MyProject');

    await waitFor(() => {
      expect(useProjectStore.getState().currProjectName).toBe('MyProject');
      expect(useProjectStore.getState().currProject).toEqual(project);
    });
  });

  it('keeps the dialog open when creation fails', async () => {
    mockPostJson.mockResolvedValueOnce(null);

    await openDialogAndSubmit('MyProject');

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(screen.queryByRole('dialog')).toBeTruthy();
  });
});
