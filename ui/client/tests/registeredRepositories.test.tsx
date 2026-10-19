import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, waitFor, cleanup, act, renderHook } from '@testing-library/react';
import { SimpleTreeView } from '@mui/x-tree-view';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockGetJson = vi.fn();
vi.mock('../src/io/entrypoint', () => ({
  getJson: (...args: any[]) => mockGetJson(...args),
  postJson: vi.fn(),
}));

vi.mock('../src/components/repo/UpdateRepositoriesButton', () => ({
  UpdateRepositoriesButton: () => null,
}));

const { useRegisteredRepositories } = await import('../src/components/repo/useRegisteredRepositories');
const { RegisteredRepositories } = await import('../src/components/repo/RegisteredRepositories');
const { idFromRepoDomId } = await import('../src/shared/idDocId');

const REPO = {
  datasourceName: 'YAVNEEL',
  resource: '/repos/meteo/YAVNEEL.json',
  dataFormat: 'JSON_dict',
  toolkit: 'heradata',
  version: [0, 0, 1],
};

afterEach(() => { cleanup(); });
beforeEach(() => { vi.clearAllMocks(); });

describe('useRegisteredRepositories', () => {
  it('asks the registered entrypoint once, with no parameters', async () => {
    mockGetJson.mockResolvedValue({ repositories: [REPO] });

    const { result } = renderHook(() => useRegisteredRepositories());

    await waitFor(() => {
      expect(result.current.repositories).toEqual([REPO]);
    });
    expect(mockGetJson).toHaveBeenCalledTimes(1);
    expect(mockGetJson.mock.calls[0][0]).toBe('/repository/registered');
    expect(mockGetJson.mock.calls[0][1]).toEqual({});
  });

  it('is empty when the call failed', async () => {
    mockGetJson.mockResolvedValue(null);

    const { result } = renderHook(() => useRegisteredRepositories());

    await waitFor(() => {
      expect(mockGetJson).toHaveBeenCalledTimes(1);
    });
    expect(result.current.repositories).toEqual([]);
  });

  it('reads again when asked', async () => {
    mockGetJson.mockResolvedValue({ repositories: [] });

    const { result } = renderHook(() => useRegisteredRepositories());
    await waitFor(() => { expect(mockGetJson).toHaveBeenCalledTimes(1); });

    mockGetJson.mockResolvedValue({ repositories: [REPO] });
    await act(async () => { await result.current.fetchRepositories(); });

    expect(result.current.repositories).toEqual([REPO]);
  });
});

describe('RegisteredRepositories', () => {
  it('a row opens by its path, not by its name', async () => {
    mockGetJson.mockResolvedValue({ repositories: [REPO] });

    render(
      <SimpleTreeView defaultExpandedItems={['registered-repos']}>
        <RegisteredRepositories />
      </SimpleTreeView>,
    );

    const label = await screen.findByText('/YAVNEEL.json');
    const row = label.closest('[role="treeitem"]') as HTMLElement;
    expect(idFromRepoDomId(row.id)).toBe(REPO.resource);
  });

  it('says so when nothing is registered', async () => {
    mockGetJson.mockResolvedValue({ repositories: [] });

    render(
      <SimpleTreeView defaultExpandedItems={['registered-repos']}>
        <RegisteredRepositories />
      </SimpleTreeView>,
    );

    expect(await screen.findByText('No registered repositories')).toBeTruthy();
  });
});
