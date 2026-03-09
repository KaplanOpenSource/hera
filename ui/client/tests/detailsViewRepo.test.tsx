import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockExecPython = vi.fn();
vi.mock('../src/io/execPython', () => ({
  execPython: (...args: any[]) => mockExecPython(...args),
}));

vi.mock('../src/components/details/RepoTreeAddButton', () => ({
  RepoTreeAddButton: () => null,
}));

const { DetailsViewRepo } = await import('../src/components/details/DetailsViewRepo');

afterEach(() => { cleanup(); });
beforeEach(() => { vi.clearAllMocks(); });

describe('DetailsViewRepo', () => {
  it('calls execPython with repo path on mount', async () => {
    mockExecPython.mockResolvedValueOnce({
      data: { json: { key1: 'value1' } },
      problem: undefined,
    });

    await act(async () => {
      render(<DetailsViewRepo repoPath="/path/to/repo.json" />);
    });

    await waitFor(() => {
      expect(mockExecPython).toHaveBeenCalledTimes(1);
      const code = mockExecPython.mock.calls[0][0];
      expect(code).toContain("open('/path/to/repo.json'");
    });
  });

  it('does not call execPython for temp repos', async () => {
    await act(async () => {
      render(<DetailsViewRepo repoPath="/path/*Temp Repository*/repo.json" />);
    });

    expect(mockExecPython).not.toHaveBeenCalled();
  });

  it('handles null json response', async () => {
    mockExecPython.mockResolvedValueOnce({
      data: { json: null },
      problem: undefined,
    });

    await act(async () => {
      render(<DetailsViewRepo repoPath="/path/to/repo.json" />);
    });

    await waitFor(() => {
      expect(mockExecPython).toHaveBeenCalledTimes(1);
    });
    // Should not crash
  });

  it('handles delayed response', async () => {
    let resolveExec!: (v: any) => void;
    mockExecPython.mockImplementationOnce(() =>
      new Promise(r => { resolveExec = r; })
    );

    await act(async () => {
      render(<DetailsViewRepo repoPath="/path/to/repo.json" />);
    });

    expect(mockExecPython).toHaveBeenCalledTimes(1);

    await act(async () => {
      resolveExec({ data: { json: { key1: 'value1' } }, problem: undefined });
    });

    // After resolution, component renders tree data without crashing
    await waitFor(() => {
      expect(screen.getByText('repo.json')).toBeTruthy();
    });
  });

  it('displays tree keys after loading', async () => {
    mockExecPython.mockResolvedValueOnce({
      data: { json: { myDataSource: { path: '/data/file.csv' } } },
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
