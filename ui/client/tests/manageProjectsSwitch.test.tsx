import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const { mockApiCall } = await import('./helpers/mockApiCall');
const { FOLDERS, openDialog, rowNames, rowOf } = await import('./helpers/manageProjects');
const { useProjectStore } = await import('../src/stores/useProjectStore');

// The switch button on each row of the manage projects dialog.

beforeEach(() => {
  vi.clearAllMocks();
  const { getJson } = mockApiCall();
  getJson.mockImplementation(async (path: string) => {
    if (path === '/project/folders') {
      return FOLDERS;
    }
    return { files: [{ name: 'data.parquet', isFolder: false }] };
  });
  useProjectStore.setState({ currProjectName: 'Alpha' });
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

const switchButton = (name: string) => {
  return within(rowOf(name)).getByTestId('manage-projects-switch') as HTMLButtonElement;
};

describe('ManageProjectsDialog switch', () => {
  it('marks the current project and turns its button off', async () => {
    await openWithRows();

    expect(switchButton('Alpha').textContent).toBe('Current');
    expect(switchButton('Alpha').disabled).toBe(true);
    expect(switchButton('Beta').textContent).toBe('Switch');
    expect(switchButton('Beta').disabled).toBe(false);
  });

  it('goes to the project and closes the dialog', async () => {
    await openWithRows();

    await act(async () => { fireEvent.click(switchButton('Beta')); });

    expect(screen.getByTestId('location').textContent).toBe('/Beta');
    await waitFor(() => {
      expect(screen.queryByText('Manage projects', { selector: 'h2' })).toBeNull();
    });
  });
});
