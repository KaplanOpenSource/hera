import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetchPython = vi.fn();
vi.mock('../src/io/fetchPython', () => ({
  fetchPython: (...args: any[]) => mockFetchPython(...args),
}));

const { AddDocumentButton } = await import('../src/components/project/AddDocumentButton');
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
    toolkits: [
      { toolkit: 'LSM', cls: 'lsm.cls' },
      { toolkit: 'GIS', cls: 'gis.cls' },
    ],
  });
});

describe('AddDocumentButton', () => {
  it('opens dialog when clicked', async () => {
    render(<AddDocumentButton />);
    const btn = screen.getByRole('button', { name: /add document/i });
    await act(async () => { fireEvent.click(btn); });
    expect(await screen.findByRole('dialog')).toBeTruthy();
  });

  it('creates a regular document with name and resource', async () => {
    const updatedProject = {
      name: 'TestProject',
      documents: [{ _id: { $oid: 'new1' }, desc: { datasourceName: 'MyDoc' }, type: 'T', _cls: 'M', resource: 'myres', dataFormat: 'string', projectName: 'TestProject' }],
    };
    mockFetchPython.mockResolvedValueOnce({ data: { project: updatedProject }, problem: undefined });

    render(<AddDocumentButton />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /add document/i }));
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
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
      const code = mockFetchPython.mock.calls[0][0].code;
      expect(code).toContain("Measurements_Collection().addDocument(");
      expect(code).toContain("'TestProject'");
      expect(code).toContain("resource='myres'");
      expect(code).toContain('"datasourceName":"MyDoc"');
    });
  });

  it('creates an agent document with effects resource', async () => {
    const updatedProject = {
      name: 'TestProject',
      documents: [{ _id: { $oid: 'ag1' }, desc: { datasourceName: 'Agent1' }, type: 'ToolkitDataSource', _cls: 'M', resource: { effects: {} }, dataFormat: 'JSON_DICT', projectName: 'TestProject' }],
    };
    mockFetchPython.mockResolvedValueOnce({ data: { project: updatedProject }, problem: undefined });

    render(<AddDocumentButton />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /add document/i }));
    });

    const dialog = await screen.findByRole('dialog');
    const nameInput = within(dialog).getByRole('textbox', { name: /^name$/i });
    fireEvent.change(nameInput, { target: { value: 'Agent1' } });

    // Select Agent kind
    fireEvent.click(within(dialog).getByRole('button', { name: 'Agent' }));

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      const code = mockFetchPython.mock.calls[0][0].code;
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
    mockFetchPython.mockResolvedValueOnce({ data: { project: updatedProject }, problem: undefined });

    render(<AddDocumentButton />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /add document/i }));
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
    mockFetchPython.mockResolvedValueOnce({ data: undefined, problem: 'exec failed' });

    render(<AddDocumentButton />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /add document/i }));
    });

    const dialog = await screen.findByRole('dialog');
    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'Fail' } });

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });

    // Store should still only have the config document
    expect(useProjectStore.getState().currProject?.documents).toEqual([configDoc]);
  });

  it('handles delayed execPython response', async () => {
    let resolveExec!: (v: any) => void;
    mockFetchPython.mockImplementationOnce(() =>
      new Promise(r => { resolveExec = r; })
    );

    render(<AddDocumentButton />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /add document/i }));
    });

    const dialog = await screen.findByRole('dialog');
    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'Delayed' } });

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      expect(mockFetchPython).toHaveBeenCalledTimes(1);
    });
    // Store not yet updated — still only has the config document
    expect(useProjectStore.getState().currProject?.documents).toEqual([configDoc]);

    const updatedProject = {
      name: 'TestProject',
      documents: [{ _id: { $oid: 'd1' }, desc: { datasourceName: 'Delayed' }, type: 'T', _cls: 'M', resource: '', dataFormat: 'string', projectName: 'TestProject' }],
    };
    await act(async () => {
      resolveExec({ data: { project: updatedProject }, problem: undefined });
    });

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.documents).toHaveLength(1);
    });
  });

  it('creates a notebook document with resource derived from name', async () => {
    const updatedProject = {
      name: 'TestProject',
      documents: [configDoc, { _id: { $oid: 'nb1' }, desc: { datasourceName: 'MyNotebook' }, type: 'notebook', _cls: 'Cache', resource: '/tmp/testproject/notebooks/MyNotebook.ipynb', dataFormat: 'JSON_dict', projectName: 'TestProject' }],
    };
    mockFetchPython.mockResolvedValueOnce({ data: { project: updatedProject }, problem: undefined });

    render(<AddDocumentButton />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /add document/i }));
    });

    const dialog = await screen.findByRole('dialog');
    fireEvent.click(within(dialog).getByRole('button', { name: 'Notebook' }));

    const nameInput = within(dialog).getByRole('textbox', { name: /^name$/i });
    fireEvent.change(nameInput, { target: { value: 'MyNotebook' } });

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      const code = mockFetchPython.mock.calls[0][0].code;
      expect(code).toContain('Cache_Collection().addDocument');
      expect(code).toContain('resource="/tmp/testproject/notebooks/MyNotebook.ipynb"');
      expect(code).toContain('type="notebook"');
      expect(code).toContain('if not notebook_path.exists()');
    });
  });

  it('shows resource field with auto-generated path when notebook is selected', async () => {
    render(<AddDocumentButton />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /add document/i }));
    });

    const dialog = await screen.findByRole('dialog');
    fireEvent.click(within(dialog).getByRole('button', { name: 'Notebook' }));

    const nameInput = within(dialog).getByRole('textbox', { name: /^name$/i });
    fireEvent.change(nameInput, { target: { value: 'TestNB' } });

    const resourceInput = within(dialog).getByRole('textbox', { name: /resource/i });
    expect(resourceInput).toHaveProperty('disabled', true);
    expect((resourceInput as HTMLInputElement).value).toBe('/tmp/testproject/notebooks/TestNB.ipynb');
  });

  it('shows helper text when notebook is selected', async () => {
    render(<AddDocumentButton />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /add document/i }));
    });

    const dialog = await screen.findByRole('dialog');
    fireEvent.click(within(dialog).getByRole('button', { name: 'Notebook' }));

    expect(within(dialog).getByText(/already exists at this path/i)).toBeTruthy();
  });

  it('hides class/toolkit when notebook kind is selected', async () => {
    render(<AddDocumentButton />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /add document/i }));
    });

    const dialog = await screen.findByRole('dialog');
    // Regular mode shows the Toolkit autocomplete
    expect(within(dialog).queryByRole('combobox', { name: /toolkit/i })).toBeTruthy();

    fireEvent.click(within(dialog).getByRole('button', { name: 'Notebook' }));

    // Notebook mode hides the Toolkit autocomplete
    expect(within(dialog).queryByRole('combobox', { name: /toolkit/i })).toBeNull();
  });

  it('includes toolkit in desc when selected', async () => {
    mockFetchPython.mockResolvedValueOnce({
      data: { project: { name: 'TestProject', documents: [] } },
      problem: undefined,
    });

    render(<AddDocumentButton toolkit={{ toolkit: 'LSM', cls: 'lsm.cls' }} />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /add document/i }));
    });

    const dialog = await screen.findByRole('dialog');
    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'WithTK' } });

    await act(async () => {
      fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
    });

    await waitFor(() => {
      const code = mockFetchPython.mock.calls[0][0].code;
      expect(code).toContain('"toolkit":"LSM"');
    });
  });
});
