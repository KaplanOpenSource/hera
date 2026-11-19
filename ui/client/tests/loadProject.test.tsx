import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { screen, fireEvent, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

import {
  LOADED,
  SAMPLE_ZIP,
  brokenZip,
  clickLoad,
  loadButton,
  pickFile,
  renderButton,
  resetLoadTest,
  shown,
  zipFile,
} from './helpers/loadProjectDialog';
import { useProjectListStore } from '../src/stores/useProjectListStore';

// Picking a zip, the load dialog, and what Load sends. The name and folder checks are in loadProjectChecks.test.tsx.

let mockPostBytesForJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  ({ postBytesForJson: mockPostBytesForJson } = resetLoadTest());
});

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

describe('LoadProjectButton', () => {
  it('opens the zip file picker first, with no dialog', async () => {
    const clicked = vi.spyOn(HTMLInputElement.prototype, 'click').mockImplementation(() => { });
    renderButton();

    await act(async () => { fireEvent.click(screen.getByRole('button', { name: /load project/i })); });

    expect(clicked).toHaveBeenCalled();
    expect((screen.getByTestId('load-project-input') as HTMLInputElement).accept).toBe('.zip');
    expect(screen.queryByRole('dialog')).toBeNull();
  });

  it('shows the dialog with the file name once a zip is picked', async () => {
    renderButton();

    await pickFile(zipFile());

    expect(screen.getByRole('dialog')).toBeTruthy();
    expect(screen.getByTestId('load-project-file-name').textContent).toBe('Loaded.zip');
    expect(mockPostBytesForJson).not.toHaveBeenCalled();
  });

  it('shows the project name and folder read from the picked zip', async () => {
    renderButton();

    await pickFile(zipFile());

    await waitFor(() => {
      expect((screen.getByTestId('load-project-name') as HTMLInputElement).value).toBe('ZipSample');
    });
    const folder = screen.getByTestId('load-project-folder') as HTMLInputElement;
    expect(folder.value).toBe('/projects/ZipSample');
    expect(folder.readOnly).toBe(false);
  });

  it('says so when the picked file is not a project export', async () => {
    vi.spyOn(console, 'error').mockImplementation(() => { });
    renderButton();

    await pickFile(brokenZip());

    expect((await screen.findByTestId('load-project-read-error')).textContent).toContain('Broken.zip');
    expect(console.error).toHaveBeenCalled();
  });

  it('keeps Load off when the picked file is not a project export', async () => {
    vi.spyOn(console, 'error').mockImplementation(() => { });
    renderButton();

    await pickFile(brokenZip());

    await screen.findByTestId('load-project-read-error');
    expect(loadButton().disabled).toBe(true);
  });

  it('opens the file picker again on Choose another zip, and shows the new file', async () => {
    renderButton();
    await pickFile(zipFile());
    const clicked = vi.spyOn(HTMLInputElement.prototype, 'click').mockImplementation(() => { });

    await act(async () => { fireEvent.click(screen.getByRole('button', { name: /choose another zip/i })); });
    expect(clicked).toHaveBeenCalled();
    await pickFile(zipFile('Other.zip'));

    expect(screen.getByTestId('load-project-file-name').textContent).toBe('Other.zip');
  });

  it('sends the file bytes on Load, updates the list and switches to the loaded project', async () => {
    mockPostBytesForJson.mockResolvedValueOnce({
      projects: [{ name: 'Alpha' }, { name: 'Loaded' }],
      project: LOADED,
    });
    renderButton();
    await pickFile(zipFile());

    await clickLoad();

    await waitFor(() => {
      expect(shown.path).toBe('/Loaded');
    });
    const [path, params, bytes, options] = mockPostBytesForJson.mock.calls[0];
    expect(path).toBe('/project/load');
    expect(params).toEqual({ projectName: 'ZipSample', filesDirectory: '/projects/ZipSample' });
    expect(Array.from(new Uint8Array(bytes))).toEqual(Array.from(SAMPLE_ZIP));
    expect(options).toEqual({ label: 'load project Loaded.zip' });
    expect(shown.state).toEqual({ project: LOADED });
    expect(useProjectListStore.getState().projectNames).toEqual([{ name: 'Alpha' }, { name: 'Loaded' }]);
    await waitFor(() => {
      expect(screen.queryByRole('dialog')).toBeNull();
    });
  });

  it('sends the edited name and folder on Load, with the name trimmed', async () => {
    mockPostBytesForJson.mockResolvedValueOnce({
      projects: [{ name: 'Alpha' }, { name: 'Renamed' }],
      project: LOADED,
    });
    renderButton();
    await pickFile(zipFile());
    await waitFor(() => {
      expect((screen.getByTestId('load-project-name') as HTMLInputElement).value).toBe('ZipSample');
    });

    fireEvent.change(screen.getByTestId('load-project-name'), { target: { value: ' Renamed ' } });
    fireEvent.change(screen.getByTestId('load-project-folder'), { target: { value: '/projects/Renamed' } });
    await clickLoad();

    await waitFor(() => {
      expect(mockPostBytesForJson).toHaveBeenCalledTimes(1);
    });
    const [path, params, bytes] = mockPostBytesForJson.mock.calls[0];
    expect(path).toBe('/project/load');
    expect(params).toEqual({ projectName: 'Renamed', filesDirectory: '/projects/Renamed' });
    expect(Array.from(new Uint8Array(bytes))).toEqual(Array.from(SAMPLE_ZIP));
  });

  it('keeps Load off and the fields read only while the zip is still being read', async () => {
    const file = zipFile();
    vi.spyOn(file, 'arrayBuffer').mockReturnValue(new Promise<ArrayBuffer>(() => { }));
    renderButton();

    await pickFile(file);

    expect(loadButton().disabled).toBe(true);
    expect((screen.getByTestId('load-project-name') as HTMLInputElement).readOnly).toBe(true);
    expect((screen.getByTestId('load-project-folder') as HTMLInputElement).readOnly).toBe(true);
  });

  it('reads the new zip after Choose another zip, dropping the old one\'s name', async () => {
    vi.spyOn(console, 'error').mockImplementation(() => { });
    vi.spyOn(HTMLInputElement.prototype, 'click').mockImplementation(() => { });
    renderButton();
    await pickFile(zipFile());
    await waitFor(() => {
      expect((screen.getByTestId('load-project-name') as HTMLInputElement).value).toBe('ZipSample');
    });

    await act(async () => { fireEvent.click(screen.getByRole('button', { name: /choose another zip/i })); });
    await pickFile(brokenZip());

    await screen.findByTestId('load-project-read-error');
    expect((screen.getByTestId('load-project-name') as HTMLInputElement).value).toBe('');
    expect((screen.getByTestId('load-project-folder') as HTMLInputElement).value).toBe('');
    expect(loadButton().disabled).toBe(true);
  });

  it('keeps the dialog open and changes nothing when the server refused the load', async () => {
    mockPostBytesForJson.mockResolvedValueOnce(null);
    renderButton();
    await pickFile(zipFile());

    await clickLoad();

    await waitFor(() => {
      expect(mockPostBytesForJson).toHaveBeenCalledTimes(1);
    });
    expect(shown.path).toBe('/Alpha');
    expect(useProjectListStore.getState().projectNames).toEqual([{ name: 'Alpha' }]);
    expect(screen.getByRole('dialog')).toBeTruthy();
  });

  it('sends nothing on Cancel', async () => {
    renderButton();
    await pickFile(zipFile());

    await act(async () => { fireEvent.click(screen.getByRole('button', { name: /cancel/i })); });

    expect(mockPostBytesForJson).not.toHaveBeenCalled();
    expect(shown.path).toBe('/Alpha');
    await waitFor(() => {
      expect(screen.queryByRole('dialog')).toBeNull();
    });
  });

  it('shows no dialog when the file picker was closed without a file', async () => {
    renderButton();

    await act(async () => {
      fireEvent.change(screen.getByTestId('load-project-input'), { target: { files: [] } });
    });

    expect(screen.queryByRole('dialog')).toBeNull();
    expect(mockPostBytesForJson).not.toHaveBeenCalled();
  });
});
