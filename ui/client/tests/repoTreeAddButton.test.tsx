import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetchPython = vi.fn();
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: (...args: any[]) => mockFetchPython(...args),
}));

const { RepoTreeAddButton } = await import('../src/components/details/RepoTreeAddButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');

afterEach(() => { cleanup(); });

beforeEach(() => {
  vi.clearAllMocks();
  useProjectStore.setState({
    currProjectName: 'TestProject',
    currProject: {
      name: 'TestProject',
      documents: [
        {
          _id: { $oid: 'cfg1' },
          type: 'TestProject__config__',
          desc: { datasourceName: 'config', filesDirectory: '/tmp/test' },
          _cls: 'M',
          resource: '',
          dataFormat: 'string',
          projectName: 'TestProject',
        },
      ],
    },
  });
});

const tree = { datasource1: { path: '/data/file1.csv' } };

const clickAddAndConfirm = async () => {
  const wrapper = screen.getByLabelText('Add repository of data sources to project');
  await act(async () => { fireEvent.click(within(wrapper).getByRole('button')); });
  const yesBtn = await screen.findByRole('button', { name: /yes/i });
  await act(async () => { fireEvent.click(yesBtn); });
};

describe('RepoTreeAddButton', () => {
  it('is disabled when tree is undefined', () => {
    render(<RepoTreeAddButton tree={undefined} />);
    const wrapper = screen.getByLabelText('Add repository of data sources to project');
    const btn = within(wrapper).getByRole('button');
    expect((btn as HTMLButtonElement).disabled).toBe(true);
  });

  it('calls execPython with repo data on confirm', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { project: { name: 'TestProject', documents: [] } },
      problem: undefined,
    });

    render(<RepoTreeAddButton tree={tree} />);
    await clickAddAndConfirm();

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
      const code = mockFetchPython.mock.calls[0][0].code;
      expect(code).toContain('loadAllDatasourcesInRepositoryJSONToProject');
      expect(code).toContain("'TestProject'");
    });
  });

  it('updates store on success', async () => {
    const updatedProject = {
      name: 'TestProject',
      documents: [
        { _id: { $oid: 'new1' }, desc: {}, type: 'T', _cls: 'M', resource: '', dataFormat: 'string', projectName: 'TestProject' },
      ],
    };
    mockFetchPython.mockResolvedValueOnce({
      data: { project: updatedProject },
      problem: undefined,
    });

    render(<RepoTreeAddButton tree={tree} />);
    await clickAddAndConfirm();

    await waitFor(() => {
      expect(useProjectStore.getState().currProject).toEqual(updatedProject);
    });
  });

  it('handles delayed response', async () => {
    let resolveExec!: (v: any) => void;
    mockFetchPython.mockImplementationOnce(() =>
      new Promise(r => { resolveExec = r; })
    );

    render(<RepoTreeAddButton tree={tree} />);
    await clickAddAndConfirm();

    expect(mockFetchPython).toHaveBeenCalledTimes(1);

    const updatedProject = { name: 'TestProject', documents: [] };
    await act(async () => {
      resolveExec({ data: { project: updatedProject }, problem: undefined });
    });

    await waitFor(() => {
      expect(useProjectStore.getState().currProject).toEqual(updatedProject);
    });
  });
});
