import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetchPython = vi.fn();
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: (...args: any[]) => mockFetchPython(...args),
}));

const { DeleteProjectButton } = await import('../src/components/header/DeleteProjectButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');

afterEach(() => {
  cleanup();
});

beforeEach(() => {
  vi.clearAllMocks();
  useProjectStore.setState({
    projectNames: [{ name: 'Alpha' }, { name: 'Beta' }],
    currProjectName: 'Alpha',
    currProject: { name: 'Alpha', documents: [] },
    toolkits: [],
  });
});

const Wrapper = () => (
  <MemoryRouter>
    <DeleteProjectButton />
  </MemoryRouter>
);

describe('DeleteProjectButton', () => {
  it('opens confirmation dialog on click', async () => {
    render(<Wrapper />);
    const wrapper = screen.getByLabelText('Delete project');
    const btn = within(wrapper).getByRole('button');
    await act(async () => { fireEvent.click(btn); });
    expect(await screen.findByText(/are you sure/i)).toBeTruthy();
  });

  it('does not delete if confirmation text does not match', async () => {
    render(<Wrapper />);
    const wrapper = screen.getByLabelText('Delete project');
    const btn = within(wrapper).getByRole('button');
    await act(async () => { fireEvent.click(btn); });

    const input = await screen.findByRole('textbox');
    fireEvent.change(input, { target: { value: 'wrong' } });

    // Confirming with a non-matching name must not trigger a delete.
    const yesBtn = screen.getByRole('button', { name: /yes/i });
    await act(async () => { fireEvent.click(yesBtn); });

    expect(mockFetchPython).not.toHaveBeenCalled();
  });

  it('calls execPython and updates store on confirmed delete', async () => {
    const remainingProject = {
      name: 'Beta',
      documents: [{ _id: { $oid: 'x' }, type: 'Beta__config__', desc: {}, _cls: 'M', resource: '', dataFormat: 'string', projectName: 'Beta' }],
    };

    mockFetchPython.mockResolvedValueOnce({
      data: {
        projectNames: [{ name: 'Beta' }],
        project: remainingProject,
      },
      problem: undefined,
    });

    render(<Wrapper />);
    const wrapper = screen.getByLabelText('Delete project');
    const btn = within(wrapper).getByRole('button');
    await act(async () => { fireEvent.click(btn); });

    const input = await screen.findByRole('textbox');
    fireEvent.change(input, { target: { value: 'Alpha' } });

    const yesBtn = screen.getByRole('button', { name: /yes/i });
    await act(async () => { fireEvent.click(yesBtn); });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });

    const code = mockFetchPython.mock.calls[0][0].code;
    expect(code).toContain("All.getDocumentsAsDict('Alpha'");
    expect(code).toContain('deleteDocumentByID');
    // Files are kept by default (checkbox unchecked).
    expect(code).toContain('deleteFiles = False');

    await waitFor(() => {
      const state = useProjectStore.getState();
      expect(state.projectNames).toEqual([{ name: 'Beta' }]);
      expect(state.currProject?.name).toBe('Beta');
    });
  });

  it('deletes files from disk when the checkbox is ticked', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { projectNames: [{ name: 'Beta' }], project: { name: 'Beta', documents: [] } },
      problem: undefined,
    });

    render(<Wrapper />);
    const wrapper = screen.getByLabelText('Delete project');
    const btn = within(wrapper).getByRole('button');
    await act(async () => { fireEvent.click(btn); });

    fireEvent.change(await screen.findByRole('textbox'), { target: { value: 'Alpha' } });
    fireEvent.click(screen.getByRole('checkbox'));

    await act(async () => { fireEvent.click(screen.getByRole('button', { name: /yes/i })); });

    await waitFor(() => {
      const code = mockFetchPython.mock.calls[0][0].code;
      expect(code).toContain('deleteFiles = True');
      // Also removes the caseConfiguration.json and drops the folder if empty.
      expect(code).toContain("'caseConfiguration.json'");
      expect(code).toContain('os.rmdir(filesDir)');
      expect(code).toContain('if not os.listdir(filesDir)');
    });
  });

  it('handles delayed execPython response', async () => {
    let resolveExec!: (v: any) => void;
    mockFetchPython.mockImplementationOnce(() =>
      new Promise(r => { resolveExec = r; })
    );

    render(<Wrapper />);
    const wrapper = screen.getByLabelText('Delete project');
    const btn = within(wrapper).getByRole('button');
    await act(async () => { fireEvent.click(btn); });

    const input = await screen.findByRole('textbox');
    fireEvent.change(input, { target: { value: 'Alpha' } });

    const yesBtn = screen.getByRole('button', { name: /yes/i });
    await act(async () => { fireEvent.click(yesBtn); });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });
    // Store not yet updated while response pending
    expect(useProjectStore.getState().currProjectName).toBe('Alpha');

    const remainingProject = {
      name: 'Beta',
      documents: [{ _id: { $oid: 'x' }, type: 'Beta__config__', desc: {}, _cls: 'M', resource: '', dataFormat: 'string', projectName: 'Beta' }],
    };
    await act(async () => {
      resolveExec({
        data: { projectNames: [{ name: 'Beta' }], project: remainingProject },
        problem: undefined,
      });
    });

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.name).toBe('Beta');
    });
  });

  it('sends delete code that sorts config document last', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { projectNames: [], project: { name: '', documents: [] } },
    });

    render(<Wrapper />);
    const wrapper = screen.getByLabelText('Delete project');
    const btn = within(wrapper).getByRole('button');
    await act(async () => { fireEvent.click(btn); });

    const input = await screen.findByRole('textbox');
    fireEvent.change(input, { target: { value: 'Alpha' } });

    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /yes/i }));
    });

    await waitFor(() => {
      const code = mockFetchPython.mock.calls[0][0].code;
      expect(code).toContain("sorted(docs, key=lambda d: d['type'] == 'Alpha__config__')");
    });
  });
});
