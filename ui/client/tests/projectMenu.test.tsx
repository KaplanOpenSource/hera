import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, cleanup, act } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

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
    const file = new File([new Uint8Array([80, 75, 3, 4])], 'Loaded.zip', { type: 'application/zip' });
    await act(async () => { fireEvent.change(input, { target: { files: [file] } }); });

    expect(mockPostBytesForJson).toHaveBeenCalledTimes(1);
    expect(screen.queryByRole('button', { name: /load project/i })).toBeNull();
  });
});
