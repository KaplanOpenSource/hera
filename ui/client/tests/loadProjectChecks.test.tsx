import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { screen, fireEvent, waitFor, cleanup } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

import { follows, isBetween, loadButton, pickFile, renderButton, resetLoadTest, zipFile } from './helpers/loadProjectDialog';

// The load dialog's checks of the name and folder against the other projects.

let mockGetJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  ({ getJson: mockGetJson } = resetLoadTest());
});

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

describe('LoadProjectDialog checks', () => {
  it('asks for the other projects\' folders when the dialog opens', async () => {
    renderButton();

    await pickFile(zipFile());

    expect(mockGetJson).toHaveBeenCalledWith('/project/folders', {}, { label: 'project folders' });
  });

  it('blocks Load when a project with the zip\'s name exists', async () => {
    mockGetJson.mockResolvedValue({ projects: [{ name: 'ZipSample', filesDirectory: '/elsewhere' }] });
    renderButton();

    await pickFile(zipFile());

    expect((await screen.findByTestId('load-project-name-taken')).textContent).toContain('ZipSample');
    expect(isBetween('load-project-name', 'load-project-name-taken', 'load-project-folder')).toBe(true);
    expect(loadButton().disabled).toBe(true);
    expect(screen.queryByTestId('load-project-shared-folder')).toBeNull();
  });

  it('shows both the name error and the folder warning when the same project still exists', async () => {
    mockGetJson.mockResolvedValue({ projects: [{ name: 'ZipSample', filesDirectory: '/projects/ZipSample' }] });
    renderButton();

    await pickFile(zipFile());

    expect((await screen.findByTestId('load-project-name-taken')).textContent).toContain('ZipSample');
    expect(screen.getByTestId('load-project-shared-folder').textContent).toContain('ZipSample');
    expect(loadButton().disabled).toBe(true);
  });

  it('reruns both checks as the name and folder are edited', async () => {
    mockGetJson.mockResolvedValue({ projects: [
      { name: 'Alpha', filesDirectory: '/projects/Alpha' },
      { name: 'ZipSample', filesDirectory: '/projects/ZipSample' },
    ] });
    renderButton();
    await pickFile(zipFile());
    await screen.findByTestId('load-project-name-taken');

    fireEvent.change(screen.getByTestId('load-project-name'), { target: { value: 'Fresh' } });
    fireEvent.change(screen.getByTestId('load-project-folder'), { target: { value: '/projects/New' } });
    expect(screen.queryByTestId('load-project-name-taken')).toBeNull();
    expect(screen.queryByTestId('load-project-shared-folder')).toBeNull();
    expect(loadButton().disabled).toBe(false);

    fireEvent.change(screen.getByTestId('load-project-name'), { target: { value: 'Alpha' } });
    fireEvent.change(screen.getByTestId('load-project-folder'), { target: { value: '/projects/Alpha' } });
    expect(screen.getByTestId('load-project-name-taken').textContent).toContain('Alpha');
    expect(screen.getByTestId('load-project-shared-folder').textContent).toContain('Alpha');
    expect(loadButton().disabled).toBe(true);
  });

  it('warns but allows Load when other projects use the zip\'s folder', async () => {
    mockGetJson.mockResolvedValue({ projects: [
      { name: 'Alpha', filesDirectory: '/projects/ZipSample' },
      { name: 'Beta', filesDirectory: '/projects/ZipSample' },
      { name: 'Gamma', filesDirectory: '/projects/Gamma' },
    ] });
    renderButton();

    await pickFile(zipFile());

    expect((await screen.findByTestId('load-project-shared-folder')).textContent).toContain('Alpha, Beta');
    expect(screen.getByTestId('load-project-shared-folder').textContent).not.toContain('Gamma');
    expect(follows('load-project-folder', 'load-project-shared-folder')).toBe(true);
    await waitFor(() => {
      expect(loadButton().disabled).toBe(false);
    });
    expect(screen.queryByTestId('load-project-name-taken')).toBeNull();
  });

  it('does not warn about projects with no folder when the folder is cleared', async () => {
    mockGetJson.mockResolvedValue({ projects: [{ name: 'NoFolder', filesDirectory: '' }] });
    renderButton();
    await pickFile(zipFile());
    await waitFor(() => {
      expect(loadButton().disabled).toBe(false);
    });

    fireEvent.change(screen.getByTestId('load-project-folder'), { target: { value: '' } });

    expect(screen.queryByTestId('load-project-shared-folder')).toBeNull();
  });

  it('still allows Load when the folders could not be fetched', async () => {
    mockGetJson.mockResolvedValue(null);
    renderButton();

    await pickFile(zipFile());

    await waitFor(() => {
      expect(loadButton().disabled).toBe(false);
    });
    expect(screen.queryByTestId('load-project-name-taken')).toBeNull();
    expect(screen.queryByTestId('load-project-shared-folder')).toBeNull();
  });

  it('blocks Load while the name is empty or only spaces', async () => {
    renderButton();
    await pickFile(zipFile());
    await waitFor(() => {
      expect(loadButton().disabled).toBe(false);
    });

    fireEvent.change(screen.getByTestId('load-project-name'), { target: { value: '' } });
    expect(loadButton().disabled).toBe(true);

    fireEvent.change(screen.getByTestId('load-project-name'), { target: { value: '   ' } });
    expect(loadButton().disabled).toBe(true);

    fireEvent.change(screen.getByTestId('load-project-name'), { target: { value: 'Fresh' } });
    expect(loadButton().disabled).toBe(false);
  });
});
