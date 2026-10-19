import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockGetJson = vi.fn();
vi.mock('../src/io/entrypoint', () => ({
  getJson: (...args: any[]) => mockGetJson(...args),
  postJson: vi.fn(),
}));

vi.mock('../src/components/details/RepoTreeAddButton', () => ({
  RepoTreeAddButton: () => null,
}));

const { DetailsViewRepo } = await import('../src/components/details/DetailsViewRepo');
const { DetailsViewMergedRepo } = await import('../src/components/details/DetailsViewMergedRepo');

afterEach(() => { cleanup(); });
beforeEach(() => {
  vi.clearAllMocks();
  localStorage.clear();
});

const REPO_PATH = '/path/to/repo.json';

const answerWith = (json: any, resolved: any = undefined) => {
  mockGetJson.mockResolvedValueOnce({ json, resolved: resolved ?? json });
};

describe('DetailsViewRepo', () => {
  it('asks the repository entrypoint for the path on mount', async () => {
    answerWith({ key1: 'value1' });

    await act(async () => {
      render(<DetailsViewRepo repoPath={REPO_PATH} />);
    });

    await waitFor(() => {
      expect(mockGetJson).toHaveBeenCalledTimes(1);
    });
    expect(mockGetJson.mock.calls[0][0]).toBe('/repository/get');
    expect(mockGetJson.mock.calls[0][1]).toEqual({ path: REPO_PATH });
  });

  it('asks by path for a registered repository too, with no branch on the name', async () => {
    answerWith({ key1: 'value1' });

    await act(async () => {
      render(<DetailsViewRepo repoPath="/repos/meteo.json" />);
    });

    await waitFor(() => {
      expect(mockGetJson.mock.calls[0][1]).toEqual({ path: '/repos/meteo.json' });
    });
  });

  it('does not call the server for temp repos', async () => {
    await act(async () => {
      render(<DetailsViewRepo repoPath="/path/*Temp Repository*/repo.json" />);
    });

    expect(mockGetJson).not.toHaveBeenCalled();
  });

  it('survives a call that failed', async () => {
    mockGetJson.mockResolvedValueOnce(null);

    await act(async () => {
      render(<DetailsViewRepo repoPath={REPO_PATH} />);
    });

    await waitFor(() => {
      expect(mockGetJson).toHaveBeenCalledTimes(1);
    });
    expect(screen.getByText('repo.json')).toBeTruthy();
  });

  it('handles a delayed answer', async () => {
    let answer!: (v: any) => void;
    mockGetJson.mockImplementationOnce(() => new Promise(r => { answer = r; }));

    await act(async () => {
      render(<DetailsViewRepo repoPath={REPO_PATH} />);
    });

    expect(mockGetJson).toHaveBeenCalledTimes(1);

    await act(async () => {
      answer({ json: { key1: 'value1' }, resolved: { key1: 'value1' } });
    });

    await waitFor(() => {
      expect(screen.getByText('repo.json')).toBeTruthy();
    });
  });

  it('shows the resolved tree, which is what the answer carries', async () => {
    answerWith(
      { myDataSource: { resource: 'file.csv' } },
      { myDataSource: { resource: '/data/file.csv' } },
    );

    await act(async () => {
      render(<DetailsViewRepo repoPath={REPO_PATH} />);
    });

    await waitFor(() => {
      expect(screen.getByText('myDataSource')).toBeTruthy();
    });
  });
});

describe('DetailsViewMergedRepo', () => {
  it('asks the folder entrypoint for the contents', async () => {
    mockGetJson.mockResolvedValueOnce({ files: [], repoJsons: {} });

    await act(async () => {
      render(<DetailsViewMergedRepo />);
    });

    await waitFor(() => {
      expect(mockGetJson).toHaveBeenCalledTimes(1);
    });
    expect(mockGetJson.mock.calls[0][0]).toBe('/repository/jsonsInFolder');
    expect(mockGetJson.mock.calls[0][1]).toEqual({
      folder: '~/hera/repositories/',
      withContents: 'true',
    });
  });

  it('merges the files it gets back', async () => {
    mockGetJson.mockResolvedValueOnce({
      files: ['/repos/a.json', '/repos/b.json'],
      repoJsons: {
        '/repos/a.json': { ToolkitA: { config: {} } },
        '/repos/b.json': { ToolkitB: { cache: {} } },
      },
    });

    await act(async () => {
      render(<DetailsViewMergedRepo />);
    });

    await waitFor(() => {
      expect(screen.getByText('ToolkitA')).toBeTruthy();
      expect(screen.getByText('ToolkitB')).toBeTruthy();
    });
  });

  it('says so when the folder holds no repository', async () => {
    mockGetJson.mockResolvedValueOnce({ files: [], repoJsons: {} });

    await act(async () => {
      render(<DetailsViewMergedRepo />);
    });

    await waitFor(() => {
      expect(screen.getByText(/No valid repository files found/)).toBeTruthy();
    });
  });
});
