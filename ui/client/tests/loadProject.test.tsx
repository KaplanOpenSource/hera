import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, cleanup, act } from '@testing-library/react';
import { MemoryRouter, Route, Routes, useLocation } from 'react-router-dom';
import { readFileSync } from 'fs';
import { dirname, join } from 'path';
import { fileURLToPath } from 'url';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

import { mockApiCall } from './helpers/mockApiCall';

const { LoadProjectButton } = await import('../src/components/header/LoadProjectButton');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');

const LOADED = { name: 'Loaded', documents: [{ _id: { $oid: 'cfg' } }] };

let mockPostBytesForJson: ReturnType<typeof vi.fn>;
let shownPath: string;
let shownState: any;

// Records where the button navigated to, and the state it sent along.
const WhereAmI = () => {
  const location = useLocation();
  shownPath = location.pathname;
  shownState = location.state;
  return null;
};

const renderButton = () => {
  render(
    <MemoryRouter initialEntries={['/Alpha']}>
      <LoadProjectButton />
      <Routes>
        <Route path="*" element={<WhereAmI />} />
      </Routes>
    </MemoryRouter>
  );
};

const pickFile = async (file: File) => {
  const input = screen.getByTestId('load-project-input') as HTMLInputElement;
  await act(async () => { fireEvent.change(input, { target: { files: [file] } }); });
};

const clickLoad = async () => {
  await act(async () => { fireEvent.click(screen.getByRole('button', { name: 'Load' })); });
};

const SAMPLE_ZIP = readFileSync(join(dirname(fileURLToPath(import.meta.url)), 'fixtures', 'ZipSample.zip'));

const zipFile = (name = 'Loaded.zip') => {
  return new File([SAMPLE_ZIP], name, { type: 'application/zip' });
};

beforeEach(() => {
  vi.clearAllMocks();
  ({ postBytesForJson: mockPostBytesForJson } = mockApiCall());
  useProjectListStore.setState({ projectNames: [{ name: 'Alpha' }] });
  shownPath = '';
  shownState = null;
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
    expect(folder.readOnly).toBe(true);
  });

  it('says so when the picked file is not a project export', async () => {
    vi.spyOn(console, 'error').mockImplementation(() => { });
    renderButton();

    await pickFile(new File([new Uint8Array([80, 75, 3, 4])], 'Broken.zip', { type: 'application/zip' }));

    expect((await screen.findByTestId('load-project-read-error')).textContent).toContain('Broken.zip');
    expect(console.error).toHaveBeenCalled();
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
      expect(shownPath).toBe('/Loaded');
    });
    const [path, bytes, options] = mockPostBytesForJson.mock.calls[0];
    expect(path).toBe('/project/load');
    expect(Array.from(new Uint8Array(bytes))).toEqual(Array.from(SAMPLE_ZIP));
    expect(options).toEqual({ label: 'load project Loaded.zip' });
    expect(shownState).toEqual({ project: LOADED });
    expect(useProjectListStore.getState().projectNames).toEqual([{ name: 'Alpha' }, { name: 'Loaded' }]);
    await waitFor(() => {
      expect(screen.queryByRole('dialog')).toBeNull();
    });
  });

  it('keeps the dialog open and changes nothing when the server refused the load', async () => {
    mockPostBytesForJson.mockResolvedValueOnce(null);
    renderButton();
    await pickFile(zipFile());

    await clickLoad();

    await waitFor(() => {
      expect(mockPostBytesForJson).toHaveBeenCalledTimes(1);
    });
    expect(shownPath).toBe('/Alpha');
    expect(useProjectListStore.getState().projectNames).toEqual([{ name: 'Alpha' }]);
    expect(screen.getByRole('dialog')).toBeTruthy();
  });

  it('sends nothing on Cancel', async () => {
    renderButton();
    await pickFile(zipFile());

    await act(async () => { fireEvent.click(screen.getByRole('button', { name: /cancel/i })); });

    expect(mockPostBytesForJson).not.toHaveBeenCalled();
    expect(shownPath).toBe('/Alpha');
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
