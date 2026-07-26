import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetchPython = vi.fn();
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: (...args: any[]) => mockFetchPython(...args),
}));

const { CopyProjectButton } = await import('../src/components/header/CopyProjectButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');

afterEach(() => {
  cleanup();
});

beforeEach(() => {
  vi.clearAllMocks();
  useProjectStore.setState({
    projectNames: [{ name: 'Alpha' }],
    currProjectName: 'Alpha',
    currProject: { name: 'Alpha', documents: [] },
    toolkits: [],
  });
});

const Wrapper = () => (
  <MemoryRouter>
    <CopyProjectButton />
  </MemoryRouter>
);

const openDialog = async () => {
  render(<Wrapper />);
  const wrapper = screen.getByLabelText('Copy project');
  const btn = within(wrapper).getByRole('button');
  await act(async () => { fireEvent.click(btn); });
  return screen.findByRole('dialog');
};

describe('CopyProjectButton', () => {
  it('prefills a default copy name', async () => {
    const dialog = await openDialog();
    const nameInput = within(dialog).getByRole('textbox', { name: /new project name/i }) as HTMLInputElement;
    expect(nameInput.value).toBe('Alpha_copy');
  });

  it('prefills the folder with a sibling of the source files directory', async () => {
    useProjectStore.setState({
      projectNames: [{ name: 'Alpha' }],
      currProjectName: 'Alpha',
      currProject: {
        name: 'Alpha',
        documents: [{
          _id: { $oid: 'cfg' },
          type: 'Alpha__config__',
          desc: { datasourceName: 'config', filesDirectory: '/data/projects/Alpha' },
          _cls: 'Cache',
          resource: '',
          dataFormat: 'string',
          projectName: 'Alpha',
        }],
      },
      toolkits: [],
    });

    const dialog = await openDialog();
    const dirInput = within(dialog).getByRole('textbox', { name: /new files directory/i }) as HTMLInputElement;
    expect(dirInput.value).toBe('/data/projects/Alpha_copy');
  });

  it('lists the distinct folders the project files come from', async () => {
    useProjectStore.setState({
      projectNames: [{ name: 'Alpha' }],
      currProjectName: 'Alpha',
      currProject: {
        name: 'Alpha',
        documents: [
          { _id: { $oid: 'cfg' }, type: 'Alpha__config__', desc: { filesDirectory: '/data/Alpha' }, _cls: 'Cache', resource: '', dataFormat: 'string', projectName: 'Alpha' },
          { _id: { $oid: 'a' }, type: 'notebook', desc: {}, _cls: 'Cache', resource: '/data/Alpha/notebooks/x.ipynb', dataFormat: 'JSON_dict', projectName: 'Alpha' },
          { _id: { $oid: 'b' }, type: 'notebook', desc: {}, _cls: 'Cache', resource: '/data/Alpha/notebooks/y.ipynb', dataFormat: 'JSON_dict', projectName: 'Alpha' },
          { _id: { $oid: 'c' }, type: 'T', desc: {}, _cls: 'Measurements', resource: '/data/Alpha/raw/z.parquet', dataFormat: 'parquet', projectName: 'Alpha' },
          { _id: { $oid: 'd' }, type: 'agent', desc: {}, _cls: 'Cache', resource: { effects: {} }, dataFormat: 'JSON_dict', projectName: 'Alpha' },
        ],
      },
      toolkits: [],
    });

    const dialog = await openDialog();
    // Distinct parent folders, deduped and sorted; inline (agent) resource excluded.
    expect(within(dialog).getByText('/data/Alpha/notebooks')).toBeTruthy();
    expect(within(dialog).getByText('/data/Alpha/raw')).toBeTruthy();
  });

  it('copies documents and files to the new project and folder', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { projectNames: [{ name: 'Alpha' }, { name: 'Beta' }] },
      problem: undefined,
    });

    const dialog = await openDialog();
    fireEvent.change(within(dialog).getByRole('textbox', { name: /new project name/i }), { target: { value: 'Beta' } });
    fireEvent.change(within(dialog).getByRole('textbox', { name: /new files directory/i }), { target: { value: '/data/beta' } });

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /copy/i }));
    });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
      const code = mockFetchPython.mock.calls[0][0].code;
      expect(code).toContain("srcName = 'Alpha'");
      expect(code).toContain("newName = 'Beta'");
      expect(code).toContain("newDir = os.path.abspath('/data/beta')");
      expect(code).toContain('createProjectDirectory(outputPath=newDir, projectName=newName)');
      expect(code).toContain('addDocumentFromDict');
      expect(code).toContain('os.path.relpath(resAbs, srcDir)');
      expect(code).toContain('if newName in getProjectList()');
    });

    await waitFor(() => {
      expect(useProjectStore.getState().projectNames).toEqual([{ name: 'Alpha' }, { name: 'Beta' }]);
    });
  });

  it('uses a default folder when none is given', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { projectNames: [{ name: 'Alpha' }, { name: 'Alpha_copy' }] },
      problem: undefined,
    });

    const dialog = await openDialog();
    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /copy/i }));
    });

    await waitFor(() => {
      const code = mockFetchPython.mock.calls[0][0].code;
      expect(code).toContain("newDir = os.path.abspath(os.path.join(os.getcwd(), 'projects', 'Alpha_copy'))");
    });
  });
});
