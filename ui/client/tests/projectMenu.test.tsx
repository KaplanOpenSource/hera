import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, cleanup, act, waitFor } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';
import { readFileSync } from 'fs';
import { dirname, join } from 'path';
import { fileURLToPath } from 'url';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

import { mockApiCall } from './helpers/mockApiCall';

const { ProjectMenuButton } = await import('../src/components/header/ProjectMenuButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');

let mockPostBytesForJson: ReturnType<typeof vi.fn>;

const renderMenu = () => {
  render(
    <MemoryRouter initialEntries={['/Alpha']}>
      <ProjectMenuButton />
    </MemoryRouter>
  );
};

const openMenu = async () => {
  await act(async () => {
    fireEvent.click(screen.getByTestId('MoreVertIcon').closest('button')!);
  });
};

beforeEach(() => {
  vi.clearAllMocks();
  ({ postBytesForJson: mockPostBytesForJson } = mockApiCall());
  useProjectStore.setState({ currProjectName: 'Alpha', currProject: { name: 'Alpha', documents: [] } });
  useProjectListStore.setState({ projectNames: [{ name: 'Alpha' }] });
});

afterEach(() => {
  cleanup();
});

describe('ProjectMenuButton', () => {
  it('opens a popover with the four project actions', async () => {
    renderMenu();
    await openMenu();

    expect(screen.getByText('Project')).toBeTruthy();
    expect(screen.getByRole('button', { name: /add project/i })).toBeTruthy();
    expect(screen.getByRole('button', { name: /load project/i })).toBeTruthy();
    expect(screen.getByRole('button', { name: /export project/i })).toBeTruthy();
    expect(screen.getByRole('button', { name: /delete project/i })).toBeTruthy();
  });

  it('closes after a project is loaded', async () => {
    mockPostBytesForJson.mockResolvedValueOnce({
      projects: [{ name: 'Alpha' }, { name: 'Loaded' }],
      project: { name: 'Loaded', documents: [] },
    });
    renderMenu();
    await openMenu();

    const input = screen.getByTestId('load-project-input') as HTMLInputElement;
    const sample = readFileSync(join(dirname(fileURLToPath(import.meta.url)), 'fixtures', 'ZipSample.zip'));
    const file = new File([sample], 'Loaded.zip', { type: 'application/zip' });
    await act(async () => { fireEvent.change(input, { target: { files: [file] } }); });
    // Load turns on once the zip is read.
    await waitFor(() => {
      expect((screen.getByRole('button', { name: 'Load' }) as HTMLButtonElement).disabled).toBe(false);
    });
    await act(async () => { fireEvent.click(screen.getByRole('button', { name: 'Load' })); });

    expect(mockPostBytesForJson).toHaveBeenCalledTimes(1);
    expect(screen.queryByRole('button', { name: /load project/i })).toBeNull();
  });
});
