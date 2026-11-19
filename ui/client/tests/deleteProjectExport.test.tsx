import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockPushError = vi.fn();
vi.mock('../src/io/snackbar', () => ({
  pushError: (...args: any[]) => mockPushError(...args),
  pushInfo: vi.fn(),
  pushRunning: vi.fn(),
  dismiss: vi.fn(),
}));

const { mockApiCall } = await import('./helpers/mockApiCall');
const { DeleteProjectButton } = await import('../src/components/header/DeleteProjectButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');

let mockPostJson: ReturnType<typeof vi.fn>;
let mockPostBytes: ReturnType<typeof vi.fn>;
let clickedDownloads: string[];

beforeEach(() => {
  vi.clearAllMocks();
  const mocks = mockApiCall();
  ({ postJson: mockPostJson, postBytes: mockPostBytes } = mocks);
  mocks.getJson.mockResolvedValue({ files: [{ name: 'data.parquet', isFolder: false }] });
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

const tickExport = () => {
  fireEvent.click(screen.getByLabelText(/export the project as a zip first/i));
};

const nameField = () => {
  return screen.getByRole('textbox') as HTMLInputElement;
};

const exportRadio = () => {
  return screen.getByLabelText(/export the project as a zip first/i) as HTMLInputElement;
};

describe('DeleteProjectButton export first', () => {
  it('saves the zip and deletes with no name typed when files are kept', async () => {
    mockPostBytes.mockResolvedValueOnce(new Uint8Array([1, 2, 3]).buffer);
    mockPostJson.mockResolvedValueOnce({ projects: [{ name: 'Beta' }], project: { name: 'Beta', documents: [] } });
    await openDialog();

    expect(nameField().disabled).toBe(false);
    tickExport();
    expect(nameField().disabled).toBe(true);
    await clickYes();

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(mockPostBytes).toHaveBeenCalledWith('/project/export',
      { projectName: 'Alpha' }, { label: 'export project Alpha' });
    expect(clickedDownloads).toEqual(['Alpha.zip']);
    expect(mockPostJson.mock.calls[0][1]).toEqual({ projectName: 'Alpha', deleteFiles: false });
  });

  it('does not delete when the export failed', async () => {
    mockPostBytes.mockResolvedValueOnce(null);
    await openDialog();

    tickExport();
    await clickYes();

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledWith('Export failed, so the project was not deleted. Try again.');
    });
    expect(mockPostJson).not.toHaveBeenCalled();
    expect(clickedDownloads).toEqual([]);
  });

  it('still needs the name when the files are deleted too', async () => {
    await openDialog();

    tickExport();
    fireEvent.click(screen.getByLabelText(/also delete the project's files/i));
    expect(nameField().disabled).toBe(false);
    expect(exportRadio().disabled).toBe(true);
    await clickYes();

    expect(mockPostBytes).not.toHaveBeenCalled();
    expect(mockPostJson).not.toHaveBeenCalled();
  });

  it('does not export when the name is typed', async () => {
    mockPostJson.mockResolvedValueOnce({ projects: [{ name: 'Beta' }], project: { name: 'Beta', documents: [] } });
    await openDialog();

    fireEvent.change(screen.getByRole('textbox'), { target: { value: 'Alpha' } });
    await clickYes();

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(mockPostBytes).not.toHaveBeenCalled();
    expect(clickedDownloads).toEqual([]);
  });
});
