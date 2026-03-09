import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockExecPython = vi.fn();
vi.mock('../src/io/execPython', () => ({
  execPython: (...args: any[]) => mockExecPython(...args),
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

    // The Yes button should be disabled when text doesn't match
    const yesBtn = screen.getByRole('button', { name: /yes/i });
    expect((yesBtn as HTMLButtonElement).disabled).toBe(true);
  });

  it('calls execPython and updates store on confirmed delete', async () => {
    const remainingProject = {
      name: 'Beta',
      documents: [{ _id: { $oid: 'x' }, type: 'Beta__config__', desc: {}, _cls: 'M', resource: '', dataFormat: 'string', projectName: 'Beta' }],
    };

    mockExecPython.mockResolvedValueOnce({
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
      expect(mockExecPython).toHaveBeenCalledTimes(1);
    });

    const code = mockExecPython.mock.calls[0][0];
    expect(code).toContain("All.getDocumentsAsDict('Alpha'");
    expect(code).toContain('deleteDocumentByID');

    await waitFor(() => {
      const state = useProjectStore.getState();
      expect(state.projectNames).toEqual([{ name: 'Beta' }]);
      expect(state.currProject?.name).toBe('Beta');
    });
  });

  it('handles delayed execPython response', async () => {
    let resolveExec!: (v: any) => void;
    mockExecPython.mockImplementationOnce(() =>
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
      expect(mockExecPython).toHaveBeenCalledTimes(1);
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
    mockExecPython.mockResolvedValueOnce({
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
      const code = mockExecPython.mock.calls[0][0];
      expect(code).toContain("sorted(docs, key=lambda d: d['type'] == 'Alpha__config__')");
    });
  });
});
