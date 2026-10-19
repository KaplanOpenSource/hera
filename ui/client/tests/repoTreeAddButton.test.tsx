import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockPostJson = vi.fn();
vi.mock('../src/io/entrypoint', () => ({
  postJson: (...args: any[]) => mockPostJson(...args),
  getJson: vi.fn(),
}));

const mockPushError = vi.fn();
const mockPushInfo = vi.fn();
vi.mock('../src/io/snackbar', () => ({
  pushError: (...args: any[]) => mockPushError(...args),
  pushInfo: (...args: any[]) => mockPushInfo(...args),
  pushRunning: vi.fn(),
  dismiss: vi.fn(),
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

  it('sends the tree to the repository entrypoint on confirm', async () => {
    mockPostJson.mockResolvedValueOnce({ project: { name: 'TestProject', documents: [] }, added: 1, problems: [] });

    render(<RepoTreeAddButton tree={tree} />);
    await clickAddAndConfirm();

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(mockPostJson.mock.calls[0][0]).toBe('/repository/loadIntoProject');
    expect(mockPostJson.mock.calls[0][1]).toEqual({
      projectName: 'TestProject',
      repositoryJSON: tree,
      baseDir: '/tmp/test',
      overwrite: true,
    });
  });

  it('sends only what the user left visible', async () => {
    mockPostJson.mockResolvedValueOnce({ project: { name: 'TestProject', documents: [] }, added: 1, problems: [] });

    render(<RepoTreeAddButton tree={{ ...tree, hidden: { path: '/x' } }} hiddenPaths={new Set(['hidden'])} />);
    await clickAddAndConfirm();

    await waitFor(() => {
      expect(mockPostJson.mock.calls[0][1].repositoryJSON).toEqual(tree);
    });
  });

  it('stores the project the answer carries, and makes no second call', async () => {
    const updatedProject = {
      name: 'TestProject',
      documents: [
        { _id: { $oid: 'new1' }, desc: {}, type: 'T', _cls: 'M', resource: '', dataFormat: 'string', projectName: 'TestProject' },
      ],
    };
    mockPostJson.mockResolvedValueOnce({ project: updatedProject, added: 1, problems: [] });

    render(<RepoTreeAddButton tree={tree} />);
    await clickAddAndConfirm();

    await waitFor(() => {
      expect(useProjectStore.getState().currProject).toEqual(updatedProject);
    });
    expect(mockPostJson).toHaveBeenCalledTimes(1);
  });

  it('handles a delayed answer', async () => {
    let answer!: (v: any) => void;
    mockPostJson.mockImplementationOnce(() => new Promise(r => { answer = r; }));

    render(<RepoTreeAddButton tree={tree} />);
    await clickAddAndConfirm();

    expect(mockPostJson).toHaveBeenCalledTimes(1);

    const updatedProject = { name: 'TestProject', documents: [] };
    await act(async () => {
      answer({ project: updatedProject, added: 1, problems: [] });
    });

    await waitFor(() => {
      expect(useProjectStore.getState().currProject).toEqual(updatedProject);
    });
  });

  it('says nothing was added, and why, when hera skipped the tree', async () => {
    mockPostJson.mockResolvedValueOnce({
      project: { name: 'TestProject', documents: [] },
      added: 0,
      problems: ['GIS_Tiles: ImportError: Cannot locate class: x.y.TilesToolkit'],
    });

    render(<RepoTreeAddButton tree={tree} />);
    await clickAddAndConfirm();

    await waitFor(() => {
      expect(mockPushError).toHaveBeenCalledTimes(1);
    });
    expect(mockPushError.mock.calls[0][0]).toContain('nothing was added');
    expect(mockPushError.mock.calls[0][0]).toContain('GIS_Tiles');
  });
});
