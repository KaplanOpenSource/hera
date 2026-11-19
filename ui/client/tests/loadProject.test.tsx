import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, cleanup, act } from '@testing-library/react';
import { MemoryRouter, Route, Routes, useLocation } from 'react-router-dom';

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

const zipFile = () => {
  return new File([new Uint8Array([80, 75, 3, 4])], 'Loaded.zip', { type: 'application/zip' });
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
  it('opens the file picker for a zip', () => {
    const clicked = vi.spyOn(HTMLInputElement.prototype, 'click').mockImplementation(() => { });
    renderButton();

    fireEvent.click(screen.getByRole('button'));

    expect(clicked).toHaveBeenCalled();
    expect((screen.getByTestId('load-project-input') as HTMLInputElement).accept).toBe('.zip');
  });

  it('sends the file bytes, updates the list and switches to the loaded project', async () => {
    mockPostBytesForJson.mockResolvedValueOnce({
      projects: [{ name: 'Alpha' }, { name: 'Loaded' }],
      project: LOADED,
    });
    renderButton();

    await pickFile(zipFile());

    await waitFor(() => {
      expect(shownPath).toBe('/Loaded');
    });
    const [path, bytes, options] = mockPostBytesForJson.mock.calls[0];
    expect(path).toBe('/project/load');
    expect(Array.from(new Uint8Array(bytes))).toEqual([80, 75, 3, 4]);
    expect(options).toEqual({ label: 'load project Loaded.zip' });
    expect(shownState).toEqual({ project: LOADED });
    expect(useProjectListStore.getState().projectNames).toEqual([{ name: 'Alpha' }, { name: 'Loaded' }]);
  });

  it('changes nothing when the server refused the load', async () => {
    mockPostBytesForJson.mockResolvedValueOnce(null);
    renderButton();

    await pickFile(zipFile());

    await waitFor(() => {
      expect(mockPostBytesForJson).toHaveBeenCalledTimes(1);
    });
    expect(shownPath).toBe('/Alpha');
    expect(useProjectListStore.getState().projectNames).toEqual([{ name: 'Alpha' }]);
  });

  it('sends nothing when no file was picked', async () => {
    renderButton();

    await act(async () => {
      fireEvent.change(screen.getByTestId('load-project-input'), { target: { files: [] } });
    });

    expect(mockPostBytesForJson).not.toHaveBeenCalled();
  });
});
