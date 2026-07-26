import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetchPython = vi.fn();
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: (...args: any[]) => mockFetchPython(...args),
}));

const { AddProjectButton } = await import('../src/components/header/AddProjectButton');

afterEach(() => {
  cleanup();
});

beforeEach(() => {
  vi.clearAllMocks();
  // The create call resolves the project's files directory; other calls (repo tree) are harmless.
  mockFetchPython.mockResolvedValue({ data: { filesDirectory: '/data/existing' }, problem: undefined });
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

const callByLabel = (needle: string) =>
  mockFetchPython.mock.calls.find(c => (c[0].label as string)?.includes(needle))?.[0];

describe('AddProjectButton notebook detection', () => {
  it('creates the project, then scans it for notebooks in a separate call', async () => {
    await openDialogAndSubmit('MyProject');

    await waitFor(() => {
      const create = callByLabel('create project MyProject');
      const detect = callByLabel('detect notebooks MyProject');
      expect(create).toBeTruthy();
      expect(detect).toBeTruthy();
      // The create call returns the resolved files directory...
      expect(create.results).toContain('filesDirectory');
      // ...which the detection call scans for notebooks.
      expect(detect.code).toContain("filesDir = os.path.expanduser('/data/existing')");
      expect(detect.code).toContain("glob.glob(os.path.join(filesDir, '**', '*.ipynb'), recursive=True)");
      expect(detect.code).toContain('type="notebook"');
    });
  });

  it('does not scan if project creation fails', async () => {
    mockFetchPython.mockResolvedValue({ data: undefined, problem: 'boom' });

    await openDialogAndSubmit('MyProject');

    await waitFor(() => {
      expect(callByLabel('create project MyProject')).toBeTruthy();
    });
    expect(callByLabel('detect notebooks MyProject')).toBeUndefined();
  });
});
