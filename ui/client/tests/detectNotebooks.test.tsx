import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetchPython = vi.fn();
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: (...args: any[]) => mockFetchPython(...args),
}));

const mockPushInfo = vi.fn();
vi.mock('../src/io/snackbar', () => ({
  pushInfo: (...args: any[]) => mockPushInfo(...args),
}));

const { DetectNotebooksButton } = await import('../src/components/project/DetectNotebooksButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');

afterEach(() => {
  cleanup();
});

const configDoc = {
  _id: { $oid: 'cfg1' },
  desc: { datasourceName: 'config', filesDirectory: '/tmp/testproject' },
  type: 'TestProject__config__',
  _cls: 'Cache',
  resource: '',
  dataFormat: 'string',
  projectName: 'TestProject',
};

beforeEach(() => {
  vi.clearAllMocks();
  useProjectStore.setState({
    projectNames: [{ name: 'TestProject' }],
    currProjectName: 'TestProject',
    currProject: { name: 'TestProject', documents: [configDoc] },
  });
});

describe('DetectNotebooksButton', () => {
  it('scans the files directory for .ipynb files and skips checkpoints', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { project: { name: 'TestProject', documents: [configDoc] }, addedCount: 0 },
      problem: undefined,
    });

    render(<DetectNotebooksButton />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /detect notebooks/i }));
    });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
      const code = mockFetchPython.mock.calls[0][0].code;
      expect(code).toContain("glob.glob(os.path.join(filesDir, '**', '*.ipynb'), recursive=True)");
      expect(code).toContain(".ipynb_checkpoints");
      expect(code).toContain('type="notebook"');
      expect(code).toContain("'TestProject'");
    });
  });

  it('updates the store and reports the count when notebooks are found', async () => {
    const newNotebook = { _id: { $oid: 'nb1' }, desc: { datasourceName: 'Found' }, type: 'notebook', _cls: 'Cache', resource: '/tmp/testproject/Found.ipynb', dataFormat: 'JSON_dict', projectName: 'TestProject' };
    mockFetchPython.mockResolvedValueOnce({
      data: { project: { name: 'TestProject', documents: [configDoc, newNotebook] }, addedCount: 1 },
      problem: undefined,
    });

    const onDetected = vi.fn();
    render(<DetectNotebooksButton onDetected={onDetected} />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /detect notebooks/i }));
    });

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.documents).toHaveLength(2);
      expect(mockPushInfo).toHaveBeenCalledWith('Detected 1 new notebook');
      expect(onDetected).toHaveBeenCalled();
    });
  });

  it('reports when no new notebooks are found', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { project: { name: 'TestProject', documents: [configDoc] }, addedCount: 0 },
      problem: undefined,
    });

    render(<DetectNotebooksButton />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /detect notebooks/i }));
    });

    await waitFor(() => {
      expect(mockPushInfo).toHaveBeenCalledWith('No new notebooks found');
    });
  });

  it('does not update the store on error', async () => {
    mockFetchPython.mockResolvedValueOnce({ data: undefined, problem: 'exec failed' });

    render(<DetectNotebooksButton />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /detect notebooks/i }));
    });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });
    expect(useProjectStore.getState().currProject?.documents).toEqual([configDoc]);
    expect(mockPushInfo).not.toHaveBeenCalled();
  });
});
