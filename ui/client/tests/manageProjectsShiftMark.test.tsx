import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { screen, fireEvent, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const { mockApiCall } = await import('./helpers/mockApiCall');
const { openDialog, rowNames } = await import('./helpers/manageProjects');

// Shift-click marks every shown row between the last clicked row and this one.

const NAMES = ['Alpha', 'Beta', 'Gamma', 'Delta', 'Echo'];

beforeEach(() => {
  vi.clearAllMocks();
  const { getJson } = mockApiCall();
  const projects = [];
  for (const name of NAMES) {
    projects.push({ name, filesDirectory: `/data/${name}`, folderExists: true });
  }
  getJson.mockImplementation(async (path: string) => {
    if (path === '/project/folders') {
      return { projects };
    }
    return { files: [{ name: 'data.parquet', isFolder: false }] };
  });
});

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

const openWithRows = async () => {
  await openDialog();
  await waitFor(() => {
    expect(rowNames()).toHaveLength(NAMES.length);
  });
};

const box = (name: string) => {
  return screen.getByLabelText(`Mark ${name}`) as HTMLInputElement;
};

const click = async (name: string, shiftKey = false) => {
  await act(async () => { fireEvent.click(box(name), { shiftKey }); });
};

const markedNames = () => {
  const names: string[] = [];
  for (const name of NAMES) {
    if (box(name).checked) {
      names.push(name);
    }
  }
  return names;
};

describe('ManageProjectsDialog shift marking', () => {
  it('marks every row between the last click and the shift-click', async () => {
    await openWithRows();

    await click('Beta');
    await click('Delta', true);

    expect(markedNames()).toEqual(['Beta', 'Gamma', 'Delta']);
  });

  it('works upwards too', async () => {
    await openWithRows();

    await click('Delta');
    await click('Alpha', true);

    expect(markedNames()).toEqual(['Alpha', 'Beta', 'Gamma', 'Delta']);
  });

  it('clears the range when the shift-clicked row was marked', async () => {
    await openWithRows();
    await click('Alpha');
    await click('Echo', true);

    await click('Beta');
    await click('Delta', true);

    expect(markedNames()).toEqual(['Alpha', 'Echo']);
  });

  it('marks only the clicked row when nothing was clicked before', async () => {
    await openWithRows();

    await click('Gamma', true);

    expect(markedNames()).toEqual(['Gamma']);
  });

  it('skips rows the search hides', async () => {
    await openWithRows();
    const search = screen.getByTestId('manage-projects-search');
    // Only Alpha and Delta hold an "l".
    await act(async () => { fireEvent.change(search, { target: { value: 'l' } }); });

    await click('Alpha');
    await click('Delta', true);

    await act(async () => { fireEvent.change(search, { target: { value: '' } }); });
    expect(markedNames()).toEqual(['Alpha', 'Delta']);
  });
});
