import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { screen, fireEvent, waitFor, cleanup, act } from '@testing-library/react';
import { unzipSync } from 'fflate';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockPushError = vi.fn();
const mockPushInfo = vi.fn();
vi.mock('../src/io/snackbar', () => ({
  pushError: (...args: any[]) => mockPushError(...args),
  pushInfo: (...args: any[]) => mockPushInfo(...args),
  pushRunning: vi.fn(),
  dismiss: vi.fn(),
}));

const { mockApiCall } = await import('./helpers/mockApiCall');
const { openDialog, rowNames } = await import('./helpers/manageProjects');
const { useProjectStore } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');

// The manage projects dialog: deleting the marked projects after saving them in projects.zip.

// Alpha has an empty folder, Beta has files, Gamma has no folder, Delta's folder is missing.
const FOLDERS = {
  projects: [
    { name: 'defaultProject', filesDirectory: '/data/default', folderState: 'hasFiles' },
    { name: 'Alpha', filesDirectory: '/data/alpha', folderState: 'empty' },
    { name: 'Beta', filesDirectory: '/data/beta', folderState: 'hasFiles' },
    { name: 'Gamma', filesDirectory: '', folderState: 'noFolder' },
    { name: 'Delta', filesDirectory: '/data/delta', folderState: 'missing' },
  ],
};

let mockPostBytes: ReturnType<typeof vi.fn>;
let mockPostJson: ReturnType<typeof vi.fn>;
let savedBlobs: Blob[];
let calls: string[];
let serverNames: string[];

beforeEach(() => {
  vi.clearAllMocks();
  const { getJson, postBytes, postJson } = mockApiCall();
  mockPostBytes = postBytes;
  mockPostJson = postJson;
  getJson.mockImplementation(async (path: string) => {
    if (path === '/project/folders') {
      return FOLDERS;
    }
    return null;
  });
  calls = [];
  mockPostBytes.mockImplementation(async (_path: string, body: { projectName: string }) => {
    calls.push(`export ${body.projectName}`);
    return new TextEncoder().encode(`zip of ${body.projectName}`).buffer;
  });
  serverNames = ['defaultProject', 'Alpha', 'Beta', 'Gamma', 'Delta'];
  mockPostJson.mockImplementation(async (_path: string, body: { projectName: string }) => {
    calls.push(`delete ${body.projectName}`);
    serverNames = serverNames.filter((name) => name !== body.projectName);
    const projects = serverNames.map((name) => ({ name }));
    return { projects, project: { name: serverNames[0], documents: [] } };
  });
  savedBlobs = [];
  URL.createObjectURL = vi.fn((blob: Blob) => {
    savedBlobs.push(blob);
    return 'blob:zip';
  }) as any;
  URL.revokeObjectURL = vi.fn() as any;
  vi.spyOn(HTMLAnchorElement.prototype, 'click').mockImplementation(function (this: HTMLAnchorElement) {
    calls.push(`save ${this.download}`);
  });
  useProjectStore.setState({ currProjectName: 'Beta', currProject: { name: 'Beta', documents: [] } });
  useProjectListStore.setState({ projectNames: serverNames.map((name) => ({ name })) });
});

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

// Opens the dialog, waits for the rows, and marks the given ones.
const openAndMark = async (names: string[]) => {
  await openDialog();
  await waitFor(() => {
    expect(rowNames()).toHaveLength(4);
  });
  for (const name of names) {
    fireEvent.click(screen.getByLabelText(`Mark ${name}`));
  }
};

const deleteButton = () => {
  return screen.getByTestId('manage-projects-delete') as HTMLButtonElement;
};

const confirmButton = () => {
  return screen.getByTestId('delete-projects-confirm') as HTMLButtonElement;
};

const openConfirm = async () => {
  await act(async () => { fireEvent.click(deleteButton()); });
};

const typeCount = (text: string) => {
  fireEvent.change(screen.getByTestId('delete-projects-count'), { target: { value: text } });
};

const confirmDelete = async (count: string) => {
  await openConfirm();
  typeCount(count);
  await act(async () => { fireEvent.click(confirmButton()); });
};

describe('ManageProjectsDialog delete confirm', () => {
  it('is off when nothing is marked, and only opens the confirm', async () => {
    await openAndMark([]);
    expect(deleteButton().disabled).toBe(true);

    fireEvent.click(screen.getByLabelText('Mark Alpha'));
    fireEvent.click(screen.getByLabelText('Mark Beta'));
    expect(deleteButton().textContent).toBe('Delete 2 projects...');
    await openConfirm();

    expect(screen.getByText('Delete 2 projects?')).toBeTruthy();
    expect(mockPostBytes).not.toHaveBeenCalled();
    expect(mockPostJson).not.toHaveBeenCalled();
  });

  it('Cancel exports and deletes nothing', async () => {
    await openAndMark(['Alpha']);
    await openConfirm();
    typeCount('1');

    await act(async () => { fireEvent.click(screen.getByRole('button', { name: 'Cancel' })); });

    expect(screen.queryByTestId('delete-projects-confirm')).toBeNull();
    expect(calls).toEqual([]);
    expect(rowNames()).toHaveLength(4);
  });

  it('lists the marked projects in their groups', async () => {
    await openAndMark(['Alpha', 'Beta', 'Gamma', 'Delta']);
    await openConfirm();

    const empty = screen.getByTestId('delete-projects-empty');
    expect(empty.textContent).toContain('their files are deleted too');
    expect(empty.textContent).toContain('Alpha');
    const withFiles = screen.getByTestId('delete-projects-with-files');
    expect(withFiles.textContent).toContain('their files are kept');
    expect(withFiles.textContent).toContain('Beta - /data/beta');
    const noFolder = screen.getByTestId('delete-projects-no-folder');
    expect(noFolder.textContent).toContain('Gamma');
    expect(noFolder.textContent).toContain('Delta');
    expect(noFolder.textContent).not.toContain('Alpha');
  });

  it('keeps the delete button off until the typed number matches', async () => {
    await openAndMark(['Alpha', 'Beta', 'Gamma']);
    await openConfirm();
    expect(confirmButton().textContent).toBe('Export and delete 3 projects');
    expect(confirmButton().disabled).toBe(true);

    typeCount('4');
    expect(confirmButton().disabled).toBe(true);
    typeCount('3');
    expect(confirmButton().disabled).toBe(false);
  });
});

describe('ManageProjectsDialog delete', () => {
  it('saves projects.zip before deleting, then deletes one by one', async () => {
    await openAndMark(['Alpha', 'Beta', 'Gamma', 'Delta']);
    await confirmDelete('4');

    await waitFor(() => {
      expect(mockPushInfo).toHaveBeenCalledWith('Deleted 4 projects.');
    });
    expect(calls).toEqual([
      'export Alpha', 'export Beta', 'export Gamma', 'export Delta',
      'save projects.zip',
      'delete Alpha', 'delete Beta', 'delete Gamma', 'delete Delta',
    ]);
    const outer = unzipSync(new Uint8Array(await savedBlobs[0].arrayBuffer()));
    expect(Object.keys(outer).sort()).toEqual(['Alpha.zip', 'Beta.zip', 'Delta.zip', 'Gamma.zip']);
    expect(mockPushError).not.toHaveBeenCalled();
  });

  it('sends delete files only for empty-folder projects', async () => {
    await openAndMark(['Alpha', 'Beta', 'Gamma', 'Delta']);
    await confirmDelete('4');

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(4);
    });
    const label = (name: string) => {
      return { label: `delete project ${name}` };
    };
    expect(mockPostJson).toHaveBeenCalledWith('/project/delete', { projectName: 'Alpha', deleteFiles: true }, label('Alpha'));
    expect(mockPostJson).toHaveBeenCalledWith('/project/delete', { projectName: 'Beta', deleteFiles: false }, label('Beta'));
    expect(mockPostJson).toHaveBeenCalledWith('/project/delete', { projectName: 'Gamma', deleteFiles: false }, label('Gamma'));
    expect(mockPostJson).toHaveBeenCalledWith('/project/delete', { projectName: 'Delta', deleteFiles: false }, label('Delta'));
  });

  it('stops at a failed export: nothing saved, nothing deleted, the project is named', async () => {
    mockPostBytes.mockImplementation(async (_path: string, body: { projectName: string }) => {
      calls.push(`export ${body.projectName}`);
      if (body.projectName === 'Beta') {
        return null;
      }
      return new Uint8Array([1]).buffer;
    });
    await openAndMark(['Alpha', 'Beta', 'Gamma']);
    await confirmDelete('3');

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledWith('Delete stopped: Beta export failed.');
    });
    expect(calls).toEqual(['export Alpha', 'export Beta']);
    expect(mockPushInfo).not.toHaveBeenCalled();
    expect(rowNames()).toHaveLength(4);
    expect(deleteButton().disabled).toBe(false);
  });

  it('removes the deleted rows and marks, and updates the project list', async () => {
    await openAndMark(['Alpha', 'Gamma']);
    await confirmDelete('2');

    await waitFor(() => {
      expect(rowNames()).toEqual(['Beta', 'Delta']);
    });
    expect(screen.getByTestId('manage-projects-marked').textContent).toBe('0 of 2 marked');
    expect(useProjectListStore.getState().projectNames).toEqual(
      [{ name: 'defaultProject' }, { name: 'Beta' }, { name: 'Delta' }]);
    expect(useProjectStore.getState().currProjectName).toBe('Beta');
  });

  it('switches away from the current project when it is deleted', async () => {
    await openAndMark(['Alpha', 'Beta']);
    await confirmDelete('2');

    await waitFor(() => {
      expect(useProjectStore.getState().currProjectName).toBe('defaultProject');
    });
    expect(useProjectStore.getState().currProject?.name).toBe('defaultProject');
  });

  it('stops at a failed delete, names it, and keeps the ones already deleted', async () => {
    mockPostJson.mockImplementation(async (_path: string, body: { projectName: string }) => {
      calls.push(`delete ${body.projectName}`);
      if (body.projectName === 'Beta') {
        return null;
      }
      serverNames = serverNames.filter((name) => name !== body.projectName);
      return { projects: serverNames.map((name) => ({ name })), project: { name: serverNames[0], documents: [] } };
    });
    await openAndMark(['Alpha', 'Beta', 'Gamma']);
    await confirmDelete('3');

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledWith('Delete stopped: Beta delete failed. Deleted 1 of 3 projects.');
    });
    expect(calls).not.toContain('delete Gamma');
    expect(mockPushInfo).not.toHaveBeenCalled();
    expect(rowNames()).toEqual(['Beta', 'Gamma', 'Delta']);
    expect(useProjectListStore.getState().projectNames).not.toContainEqual({ name: 'Alpha' });
  });

  it('keeps the confirm open and locked while the work runs', async () => {
    let answer!: (v: ArrayBuffer) => void;
    mockPostBytes.mockImplementation(() => {
      return new Promise((r) => { answer = r; });
    });
    await openAndMark(['Gamma']);

    await confirmDelete('1');
    expect(confirmButton().textContent).toBe('Deleting...');
    expect(confirmButton().disabled).toBe(true);
    expect((screen.getByText('Cancel') as HTMLButtonElement).disabled).toBe(true);
    expect(screen.getByTestId('manage-projects-progress-text').textContent).toBe('Exporting 1 of 1: Gamma');

    await act(async () => { answer(new Uint8Array([1]).buffer); });
    await waitFor(() => {
      expect(mockPushInfo).toHaveBeenCalledWith('Deleted 1 project.');
    });
    expect(screen.queryByTestId('delete-projects-confirm')).toBeNull();
  });
});
