import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { screen, fireEvent, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const { mockApiCall } = await import('./helpers/mockApiCall');
const { openDialog, rowNames } = await import('./helpers/manageProjects');

// The "Only empty projects" switch: no documents, and no files on disk.

const PROJECTS = [
  { name: 'EmptyFolder', filesDirectory: '/data/e', folderState: 'empty', documentCount: 0 },
  { name: 'NoFolder', filesDirectory: '', folderState: 'noFolder', documentCount: 0 },
  { name: 'Missing', filesDirectory: '/data/gone', folderState: 'missing', documentCount: 0 },
  { name: 'CaseConfig', filesDirectory: '/data/c', folderState: 'onlyCaseConfig', documentCount: 0 },
  { name: 'HasFiles', filesDirectory: '/data/f', folderState: 'hasFiles', documentCount: 0 },
  { name: 'NotAFolder', filesDirectory: '/data/x', folderState: 'notAFolder', documentCount: 0 },
  { name: 'Unreadable', filesDirectory: '/data/u', folderState: 'unreadable', documentCount: 0 },
  { name: 'HasDocs', filesDirectory: '/data/d', folderState: 'empty', documentCount: 4 },
];

beforeEach(() => {
  vi.clearAllMocks();
  const { getJson } = mockApiCall();
  getJson.mockImplementation(async (path: string) => {
    if (path === '/project/folders') {
      return { projects: PROJECTS };
    }
    return null;
  });
});

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

const openWithRows = async () => {
  await openDialog();
  await waitFor(() => {
    expect(rowNames()).toHaveLength(PROJECTS.length);
  });
};

const toggle = async () => {
  await act(async () => { fireEvent.click(screen.getByLabelText('Only empty projects')); });
};

describe('ManageProjectsDialog empty filter', () => {
  it('shows only projects with no documents and an empty, case-config-only, missing or unset folder', async () => {
    await openWithRows();

    await toggle();

    expect(rowNames()).toEqual(['EmptyFolder', 'NoFolder', 'Missing', 'CaseConfig']);
  });

  it('shows all projects again when turned off', async () => {
    await openWithRows();

    await toggle();
    await toggle();

    expect(rowNames()).toHaveLength(PROJECTS.length);
  });

  it('works together with the search', async () => {
    await openWithRows();
    await toggle();

    fireEvent.change(screen.getByTestId('manage-projects-search'), { target: { value: 'folder' } });

    expect(rowNames()).toEqual(['EmptyFolder', 'NoFolder']);
  });
});

describe('ManageProjectsDialog mark all with a filter', () => {
  const markAllBox = () => {
    return screen.getByLabelText('Mark all') as HTMLInputElement;
  };

  const markedLine = () => {
    return screen.getByTestId('manage-projects-marked').textContent;
  };

  it('marks only the shown rows', async () => {
    await openWithRows();
    await toggle();

    await act(async () => { fireEvent.click(markAllBox()); });

    expect(markedLine()).toBe('4 of 8 marked');
    await toggle();
    expect((screen.getByLabelText('Mark HasFiles') as HTMLInputElement).checked).toBe(false);
    expect((screen.getByLabelText('Mark HasDocs') as HTMLInputElement).checked).toBe(false);
  });

  it('drops a hidden mark when marking the shown rows', async () => {
    await openWithRows();
    fireEvent.click(screen.getByLabelText('Mark HasDocs'));
    await toggle();

    await act(async () => { fireEvent.click(markAllBox()); });

    expect(markedLine()).toBe('4 of 8 marked');
    await toggle();
    expect((screen.getByLabelText('Mark HasDocs') as HTMLInputElement).checked).toBe(false);
  });

  it('shows the header box by the shown rows only', async () => {
    await openWithRows();
    fireEvent.click(screen.getByLabelText('Mark HasDocs'));

    await toggle();

    expect(markAllBox().checked).toBe(false);
    expect(markAllBox().getAttribute('data-indeterminate')).toBe('false');
  });
});
