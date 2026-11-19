import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const { mockApiCall } = await import('./helpers/mockApiCall');
const { DeleteProjectButton } = await import('../src/components/header/DeleteProjectButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');

let mockGetJson: ReturnType<typeof vi.fn>;
let mockPostJson: ReturnType<typeof vi.fn>;
let mockPostBytes: ReturnType<typeof vi.fn>;
let clickedDownloads: string[];

const EMPTY_NOTE = /files folder is empty, so there is nothing to delete/i;

beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson: mockGetJson, postJson: mockPostJson, postBytes: mockPostBytes } = mockApiCall());
  URL.createObjectURL = vi.fn().mockReturnValue('blob:zip') as any;
  URL.revokeObjectURL = vi.fn() as any;
  clickedDownloads = [];
  vi.spyOn(HTMLAnchorElement.prototype, 'click').mockImplementation(function (this: HTMLAnchorElement) {
    clickedDownloads.push(this.download);
  });
  useProjectStore.setState({
    currProjectName: 'Alpha',
    currProject: { name: 'Alpha', documents: [] },
  });
  useProjectListStore.setState({ projectNames: [{ name: 'Alpha' }, { name: 'Beta' }] });
});

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

const openDialog = async () => {
  render(<MemoryRouter><DeleteProjectButton /></MemoryRouter>);
  const btn = within(screen.getByLabelText('Delete project')).getByRole('button');
  await act(async () => { fireEvent.click(btn); });
  await screen.findByText(/are you sure/i);
};

const clickYes = async () => {
  await act(async () => { fireEvent.click(screen.getByRole('button', { name: /yes/i })); });
};

const filesCheckbox = () => {
  return screen.getByLabelText(/also delete the project's files/i) as HTMLInputElement;
};

const exportRadio = () => {
  return screen.getByLabelText(/export the project as a zip first/i) as HTMLInputElement;
};

describe('DeleteProjectButton empty files folder', () => {
  it('asks the server for the files folder before the dialog opens', async () => {
    mockGetJson.mockResolvedValueOnce({ files: [] });
    await openDialog();

    expect(mockGetJson).toHaveBeenCalledWith('/project/files',
      { projectName: 'Alpha', path: '' }, { label: 'project files' });
  });

  it('shows the files checkbox ticked and disabled, with a note, when the folder is empty', async () => {
    mockGetJson.mockResolvedValueOnce({ files: [] });
    await openDialog();

    expect(filesCheckbox().checked).toBe(true);
    expect(filesCheckbox().disabled).toBe(true);
    expect(screen.getByText(EMPTY_NOTE)).toBeTruthy();
    expect(exportRadio().disabled).toBe(false);
  });

  it('exports and deletes with no name typed when the folder is empty', async () => {
    mockGetJson.mockResolvedValueOnce({ files: [] });
    mockPostBytes.mockResolvedValueOnce(new Uint8Array([1, 2, 3]).buffer);
    mockPostJson.mockResolvedValueOnce({ projects: [{ name: 'Beta' }], project: { name: 'Beta', documents: [] } });
    await openDialog();

    fireEvent.click(exportRadio());
    await clickYes();

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(clickedDownloads).toEqual(['Alpha.zip']);
    // The empty folder goes too, so no stale case configuration is left.
    expect(mockPostJson.mock.calls[0][1]).toEqual({ projectName: 'Alpha', deleteFiles: true });
  });

  it('counts a folder with only the case configuration as empty', async () => {
    mockGetJson.mockResolvedValueOnce({ files: [{ name: 'caseConfiguration.json', isFolder: false }] });
    await openDialog();

    expect(filesCheckbox().checked).toBe(true);
    expect(filesCheckbox().disabled).toBe(true);
    expect(screen.getByText(EMPTY_NOTE)).toBeTruthy();
  });

  it('shows no note and leaves the checkbox free when the folder has files', async () => {
    mockGetJson.mockResolvedValueOnce({ files: [{ name: 'data.parquet', isFolder: false }] });
    await openDialog();

    expect(filesCheckbox().checked).toBe(false);
    expect(filesCheckbox().disabled).toBe(false);
    expect(screen.queryByText(EMPTY_NOTE)).toBeNull();
  });

  it('works as before when the files call fails', async () => {
    mockGetJson.mockResolvedValueOnce(null);
    await openDialog();

    expect(filesCheckbox().checked).toBe(false);
    expect(filesCheckbox().disabled).toBe(false);
    expect(screen.queryByText(EMPTY_NOTE)).toBeNull();
    fireEvent.click(filesCheckbox());
    expect(exportRadio().disabled).toBe(true);
    await clickYes();

    expect(mockPostBytes).not.toHaveBeenCalled();
    expect(mockPostJson).not.toHaveBeenCalled();
  });
});
