import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { screen, fireEvent, waitFor, within, cleanup } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const { mockApiCall } = await import('./helpers/mockApiCall');
const { FOLDERS, openDialog, rowNames, rowOf } = await import('./helpers/manageProjects');

// The manage projects dialog: the list, the search, and the folder notes.

let mockGetJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson: mockGetJson } = mockApiCall());
  mockGetJson.mockImplementation(async (path: string) => {
    if (path === '/project/folders') {
      return FOLDERS;
    }
    return null;
  });
});

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

describe('ManageProjectsDialog', () => {
  it('lists all projects but defaultProject, with the count', async () => {
    await openDialog();

    await waitFor(() => {
      expect(rowNames()).toEqual(['Alpha', 'Beta', 'Gamma']);
    });
    expect(screen.getByTestId('manage-projects-count').textContent).toBe('3 projects');
    expect(mockGetJson).toHaveBeenCalledWith('/project/folders', {}, { label: 'project folders' });
  });

  it('shows the column titles in the header', async () => {
    await openDialog();

    const headers: string[] = [];
    for (const cell of screen.getAllByRole('columnheader')) {
      headers.push(cell.textContent ?? '');
    }
    expect(headers).toEqual(['', 'Name', 'Files folder', 'Documents', 'Note', '']);
  });

  it('shows each project\'s files folder', async () => {
    await openDialog();

    await waitFor(() => {
      expect(within(rowOf('Beta')).getByText('/data/beta')).toBeTruthy();
    });
  });

  it('shows how many documents each project has', async () => {
    await openDialog();

    await waitFor(() => {
      expect(within(rowOf('Gamma')).getByTestId('manage-projects-documents').textContent).toBe('12');
    });
    expect(within(rowOf('Alpha')).getByTestId('manage-projects-documents').textContent).toBe('3');
    expect(within(rowOf('Beta')).getByTestId('manage-projects-documents').textContent).toBe('0');
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

  it('does not list any project\'s files when it opens', async () => {
    await openDialog();

    await waitFor(() => {
      expect(rowNames()).toHaveLength(3);
    });
    expect(mockGetJson).toHaveBeenCalledTimes(1);
    expect(mockGetJson).not.toHaveBeenCalledWith('/project/files', expect.anything(), expect.anything());
  });

  it('notes each folder state from the folders answer', async () => {
    mockGetJson.mockImplementation(async () => {
      return { projects: [
        { name: 'NoFolder', filesDirectory: '', folderState: 'noFolder', documentCount: 0 },
        { name: 'Missing', filesDirectory: '/data/gone', folderState: 'missing', documentCount: 0 },
        { name: 'AFile', filesDirectory: '/data/file.txt', folderState: 'notAFolder', documentCount: 0 },
        { name: 'Empty', filesDirectory: '/data/e', folderState: 'empty', documentCount: 0 },
        { name: 'CaseOnly', filesDirectory: '/data/c', folderState: 'onlyCaseConfig', documentCount: 0 },
        { name: 'HasFiles', filesDirectory: '/data/f', folderState: 'hasFiles', documentCount: 0 },
        { name: 'Locked', filesDirectory: '/data/l', folderState: 'unreadable', documentCount: 0 },
      ] };
    });
    await openDialog();
    await waitFor(() => {
      expect(rowNames()).toHaveLength(7);
    });

    const noteOf = (name: string) => {
      return within(rowOf(name)).getAllByRole('cell')[4].textContent;
    };
    expect(noteOf('NoFolder')).toBe('');
    expect(within(rowOf('NoFolder')).getByText('No directory')).toBeTruthy();
    expect(noteOf('Missing')).toBe('folder missing');
    expect(noteOf('AFile')).toBe('not a folder');
    expect(noteOf('Empty')).toBe('empty folder');
    expect(noteOf('CaseOnly')).toBe('empty folder (case config only)');
    expect(noteOf('HasFiles')).toBe('');
    expect(noteOf('Locked')).toBe('cannot read folder');
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
