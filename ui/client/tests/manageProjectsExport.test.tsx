import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
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
const { FOLDERS, openDialog, rowNames, rowOf } = await import('./helpers/manageProjects');

// The manage projects dialog: marking rows and exporting the marked ones in one zip.

let mockPostBytes: ReturnType<typeof vi.fn>;
let savedBlobs: Blob[];
let clickedDownloads: string[];

beforeEach(() => {
  vi.clearAllMocks();
  const { getJson, postBytes } = mockApiCall();
  mockPostBytes = postBytes;
  getJson.mockImplementation(async (path: string) => {
    if (path === '/project/folders') {
      return FOLDERS;
    }
    return null;
  });
  mockPostBytes.mockImplementation(async (_path: string, body: { projectName: string }) => {
    return new TextEncoder().encode(`zip of ${body.projectName}`).buffer;
  });
  savedBlobs = [];
  URL.createObjectURL = vi.fn((blob: Blob) => {
    savedBlobs.push(blob);
    return 'blob:zip';
  }) as any;
  URL.revokeObjectURL = vi.fn() as any;
  clickedDownloads = [];
  vi.spyOn(HTMLAnchorElement.prototype, 'click').mockImplementation(function (this: HTMLAnchorElement) {
    clickedDownloads.push(this.download);
  });
});

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

const openWithRows = async () => {
  await openDialog();
  await waitFor(() => {
    expect(rowNames()).toHaveLength(3);
  });
};

const markBox = (name: string) => {
  return within(rowOf(name)).getByLabelText(`Mark ${name}`) as HTMLInputElement;
};

const markAllBox = () => {
  return screen.getByLabelText('Mark all') as HTMLInputElement;
};

const markedText = () => {
  return screen.getByTestId('manage-projects-marked').textContent;
};

const exportButton = () => {
  return screen.getByTestId('manage-projects-export') as HTMLButtonElement;
};

const clickExport = async () => {
  await act(async () => { fireEvent.click(exportButton()); });
};

describe('ManageProjectsDialog marking', () => {
  it('shows 0 marked and a disabled export button at first', async () => {
    await openWithRows();

    expect(markedText()).toBe('0 of 3 marked');
    expect(exportButton().disabled).toBe(true);
    expect(markAllBox().checked).toBe(false);
  });

  it('marks and unmarks a row, and the count follows', async () => {
    await openWithRows();

    fireEvent.click(markBox('Beta'));
    expect(markBox('Beta').checked).toBe(true);
    expect(markedText()).toBe('1 of 3 marked');
    expect(exportButton().disabled).toBe(false);
    expect(exportButton().textContent).toBe('Export 1 project');

    fireEvent.click(markBox('Beta'));
    expect(markBox('Beta').checked).toBe(false);
    expect(markedText()).toBe('0 of 3 marked');
    expect(exportButton().disabled).toBe(true);
  });

  it('the header box marks all, clears all, and is half set when only some are marked', async () => {
    await openWithRows();

    fireEvent.click(markBox('Alpha'));
    expect(markAllBox().getAttribute('data-indeterminate')).toBe('true');
    expect(markAllBox().checked).toBe(false);

    fireEvent.click(markAllBox());
    expect(markedText()).toBe('3 of 3 marked');
    expect(markAllBox().checked).toBe(true);
    expect(markAllBox().getAttribute('data-indeterminate')).toBe('false');
    expect(exportButton().textContent).toBe('Export 3 projects');

    fireEvent.click(markAllBox());
    expect(markedText()).toBe('0 of 3 marked');
    expect(markBox('Gamma').checked).toBe(false);
  });

  it('search keeps the marks and the count counts all projects', async () => {
    await openWithRows();
    fireEvent.click(markBox('Alpha'));
    fireEvent.click(markBox('Gamma'));

    fireEvent.change(screen.getByTestId('manage-projects-search'), { target: { value: 'be' } });
    expect(rowNames()).toEqual(['Beta']);
    expect(markedText()).toBe('2 of 3 marked');

    fireEvent.change(screen.getByTestId('manage-projects-search'), { target: { value: '' } });
    expect(markBox('Alpha').checked).toBe(true);
    expect(markBox('Beta').checked).toBe(false);
    expect(markBox('Gamma').checked).toBe(true);
  });
});

describe('ManageProjectsDialog export', () => {
  it('saves one projects.zip holding a zip per marked project', async () => {
    await openWithRows();
    fireEvent.click(markBox('Alpha'));
    fireEvent.click(markBox('Gamma'));

    await clickExport();

    await waitFor(() => {
      expect(clickedDownloads).toEqual(['projects.zip']);
    });
    expect(mockPostBytes).toHaveBeenCalledTimes(2);
    expect(mockPostBytes).toHaveBeenCalledWith('/project/export',
      { projectName: 'Alpha' }, { label: 'export project Alpha' });
    const outer = unzipSync(new Uint8Array(await savedBlobs[0].arrayBuffer()));
    expect(Object.keys(outer).sort()).toEqual(['Alpha.zip', 'Gamma.zip']);
    expect(new TextDecoder().decode(outer['Gamma.zip'])).toBe('zip of Gamma');
    expect(mockPushInfo).toHaveBeenCalledWith('Exported 2 projects.');
    expect(mockPushError).not.toHaveBeenCalled();
  });

  it('stops at a failed export, saves nothing, and names the project', async () => {
    mockPostBytes.mockImplementation(async (_path: string, body: { projectName: string }) => {
      if (body.projectName === 'Beta') {
        return null;
      }
      return new Uint8Array([1]).buffer;
    });
    await openWithRows();
    fireEvent.click(markAllBox());

    await clickExport();

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledWith('Export stopped: Beta failed.');
    });
    expect(mockPostBytes).toHaveBeenCalledTimes(2);
    expect(mockPostBytes).not.toHaveBeenCalledWith('/project/export',
      { projectName: 'Gamma' }, { label: 'export project Gamma' });
    expect(clickedDownloads).toEqual([]);
    expect(mockPushInfo).not.toHaveBeenCalled();
    expect(exportButton().disabled).toBe(false);
  });

  it('treats a thrown export as failed', async () => {
    const errorSpy = vi.spyOn(console, 'error').mockImplementation(() => {});
    mockPostBytes.mockRejectedValueOnce(new Error('boom'));
    await openWithRows();
    fireEvent.click(markBox('Alpha'));
    fireEvent.click(markBox('Beta'));

    await clickExport();

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledWith('Export stopped: Alpha failed.');
    });
    expect(mockPostBytes).toHaveBeenCalledTimes(1);
    expect(clickedDownloads).toEqual([]);
    expect(errorSpy).toHaveBeenCalled();
  });

  it('shows a progress dialog over the list while exporting', async () => {
    let answer!: (v: ArrayBuffer) => void;
    mockPostBytes.mockImplementation(() => {
      return new Promise((r) => { answer = r; });
    });
    await openWithRows();
    fireEvent.click(markBox('Alpha'));

    await clickExport();
    expect(screen.getByTestId('export-progress-dialog')).toBeTruthy();

    await act(async () => { answer(new Uint8Array([1]).buffer); });

    expect(clickedDownloads).toEqual(['projects.zip']);
    expect(screen.queryByTestId('export-progress-dialog')).toBeNull();
    expect(exportButton().textContent).toBe('Export 1 project');
  });
});

const progressText = () => {
  return screen.queryByTestId('manage-projects-progress-text')?.textContent;
};

const progressValue = () => {
  return screen.getByTestId('manage-projects-progress-bar').getAttribute('aria-valuenow');
};

describe('ManageProjectsDialog export progress', () => {
  it('shows the bar and the count while exporting, and hides it after success', async () => {
    const answers: ((v: ArrayBuffer) => void)[] = [];
    mockPostBytes.mockImplementation(() => {
      return new Promise((r) => { answers.push(r); });
    });
    await openWithRows();
    expect(screen.queryByTestId('manage-projects-progress')).toBeNull();
    fireEvent.click(markAllBox());

    await clickExport();
    expect(progressText()).toBe('Exporting 1 of 3: Alpha');
    expect(progressValue()).toBe('0');

    await act(async () => { answers[0](new Uint8Array([1]).buffer); });
    expect(progressText()).toBe('Exporting 2 of 3: Beta');
    expect(progressValue()).toBe('33');

    await act(async () => { answers[1](new Uint8Array([1]).buffer); });
    expect(progressText()).toBe('Exporting 3 of 3: Gamma');
    expect(progressValue()).toBe('67');

    await act(async () => { answers[2](new Uint8Array([1]).buffer); });
    expect(clickedDownloads).toEqual(['projects.zip']);
    expect(screen.queryByTestId('manage-projects-progress')).toBeNull();
  });

  it('hides the bar after a failed export', async () => {
    let answer!: (v: ArrayBuffer | null) => void;
    mockPostBytes.mockImplementation(() => {
      return new Promise((r) => { answer = r; });
    });
    await openWithRows();
    fireEvent.click(markBox('Beta'));
    fireEvent.click(markBox('Gamma'));

    await clickExport();
    expect(progressText()).toBe('Exporting 1 of 2: Beta');

    await act(async () => { answer(null); });
    expect(mockPushError).toHaveBeenCalledWith('Export stopped: Beta failed.');
    expect(screen.queryByTestId('manage-projects-progress')).toBeNull();
  });
});
