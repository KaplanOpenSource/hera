import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { screen, fireEvent, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));
vi.mock('../src/io/snackbar', () => ({
  pushError: vi.fn(),
  pushInfo: vi.fn(),
  pushRunning: vi.fn(),
  dismiss: vi.fn(),
}));

const { mockApiCall } = await import('./helpers/mockApiCall');
const { openDialog, rowNames, confirmNames, confirmRowText } = await import('./helpers/manageProjects');

// Deleting projects whose folder holds only the case config, cannot be read, or is a file.

const FOLDERS = {
  projects: [
    { name: 'CaseOnly', filesDirectory: '/data/case', folderState: 'onlyCaseConfig' },
    { name: 'Locked', filesDirectory: '/data/locked', folderState: 'unreadable' },
    { name: 'AFile', filesDirectory: '/data/file.txt', folderState: 'notAFolder' },
  ],
};

let mockPostJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  const { getJson, postBytes, postJson } = mockApiCall();
  mockPostJson = postJson;
  getJson.mockImplementation(async (path: string) => {
    if (path === '/project/folders') {
      return FOLDERS;
    }
    return null;
  });
  postBytes.mockImplementation(async () => {
    return new Uint8Array([1]).buffer;
  });
  mockPostJson.mockImplementation(async () => {
    return { projects: [{ name: 'defaultProject' }], project: { name: 'defaultProject', documents: [] } };
  });
  URL.createObjectURL = vi.fn(() => 'blob:zip') as any;
  URL.revokeObjectURL = vi.fn() as any;
  vi.spyOn(HTMLAnchorElement.prototype, 'click').mockImplementation(() => {});
});

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

const openConfirmForAll = async () => {
  await openDialog();
  await waitFor(() => {
    expect(rowNames()).toHaveLength(3);
  });
  for (const name of ['CaseOnly', 'Locked', 'AFile']) {
    fireEvent.click(screen.getByLabelText(`Mark ${name}`));
  }
  await act(async () => { fireEvent.click(screen.getByTestId('manage-projects-delete')); });
};

describe('ManageProjectsDialog delete by folder state', () => {
  it('removes a case-config-only folder as empty, and keeps the others', async () => {
    await openConfirmForAll();

    expect(confirmNames()).toEqual(['CaseOnly', 'Locked', 'AFile']);
    expect(confirmRowText('CaseOnly')).toContain('Empty, removed');
    expect(confirmRowText('Locked')).toContain('/data/locked');
    expect(confirmRowText('Locked')).toContain('Kept');
    expect(confirmRowText('Locked')).toContain('cannot read folder');
    expect(confirmRowText('AFile')).toContain('Kept');
    expect(confirmRowText('AFile')).toContain('not a folder');
  });

  it('sends delete files for a case-config-only folder, but not for an unreadable one or a file', async () => {
    await openConfirmForAll();
    fireEvent.change(screen.getByTestId('delete-projects-count'), { target: { value: '3' } });
    await act(async () => { fireEvent.click(screen.getByTestId('delete-projects-confirm')); });

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(3);
    });
    const label = (name: string) => {
      return { label: `delete project ${name}` };
    };
    expect(mockPostJson).toHaveBeenCalledWith('/project/delete', { projectName: 'CaseOnly', deleteFiles: true }, label('CaseOnly'));
    expect(mockPostJson).toHaveBeenCalledWith('/project/delete', { projectName: 'Locked', deleteFiles: false }, label('Locked'));
    expect(mockPostJson).toHaveBeenCalledWith('/project/delete', { projectName: 'AFile', deleteFiles: false }, label('AFile'));
  });
});
