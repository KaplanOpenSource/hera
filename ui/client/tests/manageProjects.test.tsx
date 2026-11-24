import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const { mockApiCall } = await import('./helpers/mockApiCall');
const { ManageProjectsButton } = await import('../src/components/header/ManageProjectsButton');

// The manage projects dialog: the list, the search, and the empty-folder notes.

let mockGetJson: ReturnType<typeof vi.fn>;
let filesAnswers: { [name: string]: any };

const FOLDERS = {
  projects: [
    { name: 'defaultProject', filesDirectory: '/data/default', folderExists: true },
    { name: 'Alpha', filesDirectory: '/data/alpha', folderExists: true },
    { name: 'Beta', filesDirectory: '/data/beta', folderExists: true },
    { name: 'Gamma', filesDirectory: '/data/gamma', folderExists: true },
  ],
};

beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson: mockGetJson } = mockApiCall());
  filesAnswers = {
    Alpha: { files: [{ name: 'data.parquet', isFolder: false }] },
    Beta: { files: [{ name: 'data.parquet', isFolder: false }] },
    Gamma: { files: [{ name: 'data.parquet', isFolder: false }] },
  };
  mockGetJson.mockImplementation(async (path: string, params: any) => {
    if (path === '/project/folders') {
      return FOLDERS;
    }
    return filesAnswers[params.projectName];
  });
});

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

const openDialog = async () => {
  render(<ManageProjectsButton />);
  const btn = within(screen.getByLabelText('Manage projects')).getByRole('button');
  await act(async () => { fireEvent.click(btn); });
  await screen.findByText('Manage projects', { selector: 'h2' });
};

const rowNames = () => {
  const names: string[] = [];
  for (const row of screen.queryAllByTestId('manage-projects-row')) {
    names.push(row.querySelector('[data-testid="manage-projects-name"]')?.textContent ?? '');
  }
  return names;
};

const rowOf = (name: string) => {
  for (const row of screen.getAllByTestId('manage-projects-row')) {
    if (row.querySelector('[data-testid="manage-projects-name"]')?.textContent === name) {
      return row;
    }
  }
  throw new Error(`no row for ${name}`);
};

describe('ManageProjectsDialog', () => {
  it('lists all projects but defaultProject, with the count', async () => {
    await openDialog();

    await waitFor(() => {
      expect(rowNames()).toEqual(['Alpha', 'Beta', 'Gamma']);
    });
    expect(screen.getByTestId('manage-projects-count').textContent).toBe('3 projects');
    expect(mockGetJson).toHaveBeenCalledWith('/project/folders', {}, { label: 'project folders' });
  });

  it('shows each project\'s files folder', async () => {
    await openDialog();

    await waitFor(() => {
      expect(within(rowOf('Beta')).getByText('/data/beta')).toBeTruthy();
    });
  });

  it('filters the list by the search text', async () => {
    await openDialog();
    await waitFor(() => {
      expect(rowNames()).toHaveLength(3);
    });

    fireEvent.change(screen.getByTestId('manage-projects-search'), { target: { value: 'be' } });

    expect(rowNames()).toEqual(['Beta']);
    expect(screen.getByTestId('manage-projects-count').textContent).toBe('3 projects');
  });

  it('checks each project\'s files folder when it opens', async () => {
    await openDialog();

    await waitFor(() => {
      expect(mockGetJson).toHaveBeenCalledWith('/project/files',
        { projectName: 'Gamma', path: '' }, { label: 'project files' });
    });
    expect(mockGetJson).not.toHaveBeenCalledWith('/project/files',
      { projectName: 'defaultProject', path: '' }, { label: 'project files' });
  });

  it('notes the empty folders only', async () => {
    filesAnswers.Alpha = { files: [] };
    filesAnswers.Beta = { files: [{ name: 'caseConfiguration.json', isFolder: false }] };
    await openDialog();

    await waitFor(() => {
      expect(within(rowOf('Alpha')).queryByTestId('manage-projects-empty')).toBeTruthy();
    });
    expect(within(rowOf('Beta')).queryByTestId('manage-projects-empty')).toBeTruthy();
    expect(within(rowOf('Gamma')).queryByTestId('manage-projects-empty')).toBeNull();
    expect(screen.queryAllByTestId('manage-projects-checking')).toHaveLength(0);
  });

  it('shows checking until the folder check answers', async () => {
    let answerAlpha!: (v: any) => void;
    mockGetJson.mockImplementation(async (path: string, params: any) => {
      if (path === '/project/folders') {
        return FOLDERS;
      }
      if (params.projectName === 'Alpha') {
        return new Promise((r) => { answerAlpha = r; });
      }
      return filesAnswers[params.projectName];
    });
    await openDialog();

    await waitFor(() => {
      expect(within(rowOf('Alpha')).queryByTestId('manage-projects-checking')).toBeTruthy();
    });
    expect(within(rowOf('Beta')).queryByTestId('manage-projects-checking')).toBeNull();

    await act(async () => { answerAlpha({ files: [] }); });

    expect(within(rowOf('Alpha')).queryByTestId('manage-projects-checking')).toBeNull();
    expect(within(rowOf('Alpha')).queryByTestId('manage-projects-empty')).toBeTruthy();
  });

  it('shows no note when the folder check fails', async () => {
    filesAnswers.Alpha = null;
    await openDialog();

    await waitFor(() => {
      expect(within(rowOf('Beta')).queryByTestId('manage-projects-checking')).toBeNull();
    });
    expect(within(rowOf('Alpha')).queryByTestId('manage-projects-checking')).toBeNull();
    expect(within(rowOf('Alpha')).queryByTestId('manage-projects-empty')).toBeNull();
  });

  it('shows no note when the folder check throws', async () => {
    const errorSpy = vi.spyOn(console, 'error').mockImplementation(() => {});
    mockGetJson.mockImplementation(async (path: string, params: any) => {
      if (path === '/project/folders') {
        return FOLDERS;
      }
      if (params.projectName === 'Alpha') {
        throw new Error('boom');
      }
      return filesAnswers[params.projectName];
    });
    await openDialog();

    await waitFor(() => {
      expect(within(rowOf('Alpha')).queryByTestId('manage-projects-checking')).toBeNull();
    });
    expect(within(rowOf('Alpha')).queryByTestId('manage-projects-empty')).toBeNull();
    expect(errorSpy).toHaveBeenCalled();
  });

  it('skips the check when no folder is set or the folder is missing', async () => {
    mockGetJson.mockImplementation(async (path: string, params: any) => {
      if (path === '/project/folders') {
        return { projects: [
          { name: 'Alpha', filesDirectory: '', folderExists: false },
          { name: 'Beta', filesDirectory: '/data/gone', folderExists: false },
          { name: 'Gamma', filesDirectory: '/data/gamma', folderExists: true },
        ] };
      }
      return filesAnswers[params.projectName];
    });
    await openDialog();

    await waitFor(() => {
      expect(within(rowOf('Beta')).queryByTestId('manage-projects-missing')).toBeTruthy();
    });
    expect(within(rowOf('Alpha')).getByText('No directory')).toBeTruthy();
    expect(within(rowOf('Alpha')).queryByTestId('manage-projects-missing')).toBeNull();
    expect(within(rowOf('Alpha')).queryByTestId('manage-projects-empty')).toBeNull();
    expect(within(rowOf('Beta')).queryByTestId('manage-projects-empty')).toBeNull();
    expect(screen.queryAllByTestId('manage-projects-checking')).toHaveLength(0);
    expect(mockGetJson).toHaveBeenCalledWith('/project/files',
      { projectName: 'Gamma', path: '' }, { label: 'project files' });
    expect(mockGetJson).not.toHaveBeenCalledWith('/project/files',
      { projectName: 'Alpha', path: '' }, { label: 'project files' });
    expect(mockGetJson).not.toHaveBeenCalledWith('/project/files',
      { projectName: 'Beta', path: '' }, { label: 'project files' });
  });

  it('shows an empty list when the folders call fails', async () => {
    mockGetJson.mockResolvedValue(null);
    await openDialog();

    await waitFor(() => {
      expect(screen.getByTestId('manage-projects-count').textContent).toBe('0 projects');
    });
    expect(rowNames()).toEqual([]);
  });
});
