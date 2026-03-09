import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetchPython = vi.fn();
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: (...args: any[]) => mockFetchPython(...args),
}));

vi.mock('../src/components/details/RepoTreeAddButton', () => ({
  RepoTreeAddButton: () => null,
}));

const { DetailsViewRepo } = await import('../src/components/details/DetailsViewRepo');

afterEach(() => { cleanup(); });
beforeEach(() => { vi.clearAllMocks(); });

describe('DetailsViewRepo', () => {
  it('calls execPython with repo path on mount', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { jsonData: { key1: 'value1' } },
      problem: undefined,
    });

    await act(async () => {
      render(<DetailsViewRepo repoPath="/path/to/repo.json" />);
    });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
      const code = mockFetchPython.mock.calls[0][0].code;
      expect(code).toContain("open('/path/to/repo.json'");
    });
  });

  it('does not call execPython for temp repos', async () => {
    await act(async () => {
      render(<DetailsViewRepo repoPath="/path/*Temp Repository*/repo.json" />);
    });

    expect(mockFetchPython).not.toHaveBeenCalled();
  });

  it('handles null json response', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { jsonData: null },
      problem: undefined,
    });

    await act(async () => {
      render(<DetailsViewRepo repoPath="/path/to/repo.json" />);
    });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });
    // Should not crash
  });

  it('handles delayed response', async () => {
    let resolveExec!: (v: any) => void;
    mockFetchPython.mockImplementationOnce(() =>
      new Promise(r => { resolveExec = r; })
    );

    await act(async () => {
      render(<DetailsViewRepo repoPath="/path/to/repo.json" />);
    });

    expect(mockFetchPython).toHaveBeenCalledTimes(1);

    await act(async () => {
      resolveExec({ data: { jsonData: { key1: 'value1' } }, problem: undefined });
    });

    // After resolution, component renders tree data without crashing
    await waitFor(() => {
      expect(screen.getByText('repo.json')).toBeTruthy();
    });
  });

  it('displays tree keys after loading', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { jsonData: { myDataSource: { path: '/data/file.csv' } } },
      problem: undefined,
    });

    await act(async () => {
      render(<DetailsViewRepo repoPath="/path/to/repo.json" />);
    });

    await waitFor(() => {
      expect(screen.getByText('myDataSource')).toBeTruthy();
    });
  });
});
