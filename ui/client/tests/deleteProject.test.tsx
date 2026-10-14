import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockPostJson = vi.fn();
vi.mock('../src/io/entrypoint', () => ({
  postJson: (...args: any[]) => mockPostJson(...args),
  getJson: vi.fn(),
}));

const { DeleteProjectButton } = await import('../src/components/header/DeleteProjectButton');
const { useProjectStore, NO_PROJECT } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');

afterEach(() => {
  cleanup();
});

beforeEach(() => {
  vi.clearAllMocks();
  useProjectStore.setState({
    currProjectName: 'Alpha',
    currProject: { name: 'Alpha', documents: [] },
  });
  useProjectListStore.setState({ projectNames: [{ name: 'Alpha' }, { name: 'Beta' }] });
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

    expect(mockPostJson).not.toHaveBeenCalled();
  });

  it('calls project.delete and updates store on confirmed delete', async () => {
    const remainingProject = {
      name: 'Beta',
      documents: [{ _id: { $oid: 'x' }, type: 'Beta__config__', desc: {}, _cls: 'M', resource: '', dataFormat: 'string', projectName: 'Beta' }],
    };

    mockPostJson.mockResolvedValueOnce({
      projects: [{ name: 'Beta' }],
      project: remainingProject,
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
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });

    const [path, body] = mockPostJson.mock.calls[0];
    expect(path).toBe('/project/delete');
    // Files are kept by default (checkbox unchecked).
    expect(body).toEqual({ projectName: 'Alpha', deleteFiles: false });

    await waitFor(() => {
      const state = useProjectStore.getState();
      expect(useProjectListStore.getState().projectNames).toEqual([{ name: 'Beta' }]);
      expect(state.currProject?.name).toBe('Beta');
    });
  });

  it('deletes files from disk when the checkbox is ticked', async () => {
    mockPostJson.mockResolvedValueOnce({
      projects: [{ name: 'Beta' }],
      project: { name: 'Beta', documents: [] },
    });

    render(<Wrapper />);
    const wrapper = screen.getByLabelText('Delete project');
    const btn = within(wrapper).getByRole('button');
    await act(async () => { fireEvent.click(btn); });

    fireEvent.change(await screen.findByRole('textbox'), { target: { value: 'Alpha' } });
    fireEvent.click(screen.getByRole('checkbox'));

    await act(async () => { fireEvent.click(screen.getByRole('button', { name: /yes/i })); });

    await waitFor(() => {
      expect(mockPostJson.mock.calls[0][1].deleteFiles).toBe(true);
    });
  });

  it('handles a delayed answer', async () => {
    let resolveExec!: (v: any) => void;
    mockPostJson.mockImplementationOnce(() =>
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
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    // Store not yet updated while response pending
    expect(useProjectStore.getState().currProjectName).toBe('Alpha');

    const remainingProject = {
      name: 'Beta',
      documents: [{ _id: { $oid: 'x' }, type: 'Beta__config__', desc: {}, _cls: 'M', resource: '', dataFormat: 'string', projectName: 'Beta' }],
    };
    await act(async () => {
      resolveExec({ projects: [{ name: 'Beta' }], project: remainingProject });
    });

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.name).toBe('Beta');
    });
  });

  it('shows no project when the last one is deleted', async () => {
    mockPostJson.mockResolvedValueOnce({ projects: [], project: null });

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
      expect(useProjectListStore.getState().projectNames).toEqual([]);
      expect(useProjectStore.getState().currProject).toBeNull();
      expect(useProjectStore.getState().currProjectName).toBe(NO_PROJECT);
    });
  });

  it('does not touch the store when the delete failed', async () => {
    mockPostJson.mockResolvedValueOnce(null);

    render(<Wrapper />);
    const wrapper = screen.getByLabelText('Delete project');
    const btn = within(wrapper).getByRole('button');
    await act(async () => { fireEvent.click(btn); });

    fireEvent.change(await screen.findByRole('textbox'), { target: { value: 'Alpha' } });
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /yes/i }));
    });

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(useProjectStore.getState().currProjectName).toBe('Alpha');
    expect(useProjectListStore.getState().projectNames).toEqual([{ name: 'Alpha' }, { name: 'Beta' }]);
  });
});
