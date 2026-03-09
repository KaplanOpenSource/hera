import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockExecPython = vi.fn();
vi.mock('../src/io/execPython', () => ({
  execPython: (...args: any[]) => mockExecPython(...args),
}));

const { AddDocumentButton } = await import('../src/components/project/AddDocumentButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');

afterEach(() => {
  cleanup();
});

beforeEach(() => {
  vi.clearAllMocks();
  useProjectStore.setState({
    projectNames: [{ name: 'TestProject' }],
    currProjectName: 'TestProject',
    currProject: { name: 'TestProject', documents: [] },
    toolkits: [
      { toolkit: 'LSM', cls: 'lsm.cls' },
      { toolkit: 'GIS', cls: 'gis.cls' },
    ],
  });
});

describe('AddDocumentButton', () => {
  it('opens dialog when clicked', async () => {
    render(<AddDocumentButton />);
    const wrapper = screen.getByLabelText('Add Document');
    const btn = within(wrapper).getByRole('button');
    await act(async () => { fireEvent.click(btn); });
    expect(await screen.findByRole('dialog')).toBeTruthy();
  });

  it('creates a regular document with name and resource', async () => {
    const updatedProject = {
      name: 'TestProject',
      documents: [{ _id: { $oid: 'new1' }, desc: { datasourceName: 'MyDoc' }, type: 'T', _cls: 'M', resource: 'myres', dataFormat: 'string', projectName: 'TestProject' }],
    };
    mockExecPython.mockResolvedValueOnce({ data: updatedProject, problem: undefined });

    render(<AddDocumentButton />);
    await act(async () => {
      const w = screen.getByLabelText('Add Document');
      fireEvent.click(within(w).getByRole('button'));
    });

    const dialog = await screen.findByRole('dialog');
    const nameInput = within(dialog).getByRole('textbox', { name: /^name$/i });
    const resourceInput = within(dialog).getByRole('textbox', { name: /resource/i });

    fireEvent.change(nameInput, { target: { value: 'MyDoc' } });
    fireEvent.change(resourceInput, { target: { value: 'myres' } });

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      expect(mockExecPython).toHaveBeenCalledTimes(1);
      const code = mockExecPython.mock.calls[0][0];
      expect(code).toContain("All.addDocument('TestProject'");
      expect(code).toContain("resource='myres'");
      expect(code).toContain('"datasourceName":"MyDoc"');
    });
  });

  it('creates an agent document with effects resource', async () => {
    const updatedProject = {
      name: 'TestProject',
      documents: [{ _id: { $oid: 'ag1' }, desc: { datasourceName: 'Agent1' }, type: 'ToolkitDataSource', _cls: 'M', resource: { effects: {} }, dataFormat: 'JSON_DICT', projectName: 'TestProject' }],
    };
    mockExecPython.mockResolvedValueOnce({ data: updatedProject, problem: undefined });

    render(<AddDocumentButton />);
    await act(async () => {
      const w = screen.getByLabelText('Add Document');
      fireEvent.click(within(w).getByRole('button'));
    });

    const dialog = await screen.findByRole('dialog');
    const nameInput = within(dialog).getByRole('textbox', { name: /^name$/i });
    fireEvent.change(nameInput, { target: { value: 'Agent1' } });

    // Check the Agent checkbox
    const agentCheckbox = within(dialog).getByRole('checkbox', { name: /agent/i });
    fireEvent.click(agentCheckbox);

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      const code = mockExecPython.mock.calls[0][0];
      expect(code).toContain('resource={"effects": {}}');
      expect(code).toContain("type='ToolkitDataSource'");
    });
  });

  it('updates the store after successful creation', async () => {
    const updatedProject = {
      name: 'TestProject',
      documents: [
        { _id: { $oid: 'n1' }, desc: { datasourceName: 'NewDoc' }, type: 'T', _cls: 'M', resource: '', dataFormat: 'string', projectName: 'TestProject' },
      ],
    };
    mockExecPython.mockResolvedValueOnce({ data: updatedProject, problem: undefined });

    render(<AddDocumentButton />);
    await act(async () => {
      const w = screen.getByLabelText('Add Document');
      fireEvent.click(within(w).getByRole('button'));
    });

    const dialog = await screen.findByRole('dialog');
    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'NewDoc' } });

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      const project = useProjectStore.getState().currProject;
      expect(project?.documents).toHaveLength(1);
      expect(project?.documents[0].desc.datasourceName).toBe('NewDoc');
    });
  });

  it('does not update store on error', async () => {
    mockExecPython.mockResolvedValueOnce({ data: undefined, problem: 'exec failed' });

    render(<AddDocumentButton />);
    await act(async () => {
      const w = screen.getByLabelText('Add Document');
      fireEvent.click(within(w).getByRole('button'));
    });

    const dialog = await screen.findByRole('dialog');
    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'Fail' } });

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      expect(mockExecPython).toHaveBeenCalledTimes(1);
    });

    // Store should still have empty documents
    expect(useProjectStore.getState().currProject?.documents).toEqual([]);
  });

  it('handles delayed execPython response', async () => {
    let resolveExec!: (v: any) => void;
    mockExecPython.mockImplementationOnce(() =>
      new Promise(r => { resolveExec = r; })
    );

    render(<AddDocumentButton />);
    await act(async () => {
      const w = screen.getByLabelText('Add Document');
      fireEvent.click(within(w).getByRole('button'));
    });

    const dialog = await screen.findByRole('dialog');
    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'Delayed' } });

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      expect(mockExecPython).toHaveBeenCalledTimes(1);
    });
    // Store not yet updated
    expect(useProjectStore.getState().currProject?.documents).toEqual([]);

    const updatedProject = {
      name: 'TestProject',
      documents: [{ _id: { $oid: 'd1' }, desc: { datasourceName: 'Delayed' }, type: 'T', _cls: 'M', resource: '', dataFormat: 'string', projectName: 'TestProject' }],
    };
    await act(async () => {
      resolveExec({ data: updatedProject, problem: undefined });
    });

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.documents).toHaveLength(1);
    });
  });

  it('includes toolkit in desc when selected', async () => {
    mockExecPython.mockResolvedValueOnce({
      data: { name: 'TestProject', documents: [] },
      problem: undefined,
    });

    render(<AddDocumentButton toolkit={{ toolkit: 'LSM', cls: 'lsm.cls' }} />);
    await act(async () => {
      const w = screen.getByLabelText('Add Document');
      fireEvent.click(within(w).getByRole('button'));
    });

    const dialog = await screen.findByRole('dialog');
    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'WithTK' } });

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      const code = mockExecPython.mock.calls[0][0];
      expect(code).toContain('"toolkit":"LSM"');
    });
  });
});
