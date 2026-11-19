import { expect, vi } from 'vitest';
import { render, screen, fireEvent, waitFor, act } from '@testing-library/react';
import { MemoryRouter, Route, Routes, useLocation } from 'react-router-dom';
import { readFileSync } from 'fs';
import { dirname, join } from 'path';
import { fileURLToPath } from 'url';
import { LoadProjectButton } from '../../src/components/header/LoadProjectButton';
import { useProjectListStore } from '../../src/stores/useProjectListStore';
import { mockApiCall } from './mockApiCall';

// Shared by the load dialog tests: render the button, pick a zip, find and press Load.

export const LOADED = { name: 'Loaded', documents: [{ _id: { $oid: 'cfg' } }] };

// Where the button navigated to, and the state it sent along.
export const shown: { path: string, state: any } = { path: '', state: null };

const WhereAmI = () => {
  const location = useLocation();
  shown.path = location.pathname;
  shown.state = location.state;
  return null;
};

export const renderButton = () => {
  render(
    <MemoryRouter initialEntries={['/Alpha']}>
      <LoadProjectButton />
      <Routes>
        <Route path="*" element={<WhereAmI />} />
      </Routes>
    </MemoryRouter>
  );
};

export const pickFile = async (file: File) => {
  const input = screen.getByTestId('load-project-input') as HTMLInputElement;
  await act(async () => { fireEvent.change(input, { target: { files: [file] } }); });
};

// True when the element with id `later` comes after the one with id `earlier` in the page.
export const follows = (earlier: string, later: string) => {
  const position = screen.getByTestId(earlier).compareDocumentPosition(screen.getByTestId(later));
  return (position & Node.DOCUMENT_POSITION_FOLLOWING) !== 0;
};

export const isBetween = (before: string, middle: string, after: string) => {
  return follows(before, middle) && follows(middle, after);
};

export const loadButton = () => {
  return screen.getByRole('button', { name: 'Load' }) as HTMLButtonElement;
};

// Load turns on once the zip is read and checked.
export const clickLoad = async () => {
  await waitFor(() => {
    expect(loadButton().disabled).toBe(false);
  });
  await act(async () => { fireEvent.click(loadButton()); });
};

export const SAMPLE_ZIP = readFileSync(join(dirname(fileURLToPath(import.meta.url)), '..', 'fixtures', 'ZipSample.zip'));

export const zipFile = (name = 'Loaded.zip') => {
  return new File([SAMPLE_ZIP], name, { type: 'application/zip' });
};

export const brokenZip = (name = 'Broken.zip') => {
  return new File([new Uint8Array([80, 75, 3, 4])], name, { type: 'application/zip' });
};

// Fresh mocks before each test: one other project, Alpha, with its own folder.
export const resetLoadTest = () => {
  vi.clearAllMocks();
  const mocks = mockApiCall();
  mocks.getJson.mockResolvedValue({ projects: [{ name: 'Alpha', filesDirectory: '/projects/Alpha' }] });
  useProjectListStore.setState({ projectNames: [{ name: 'Alpha' }] });
  shown.path = '';
  shown.state = null;
  return mocks;
};
