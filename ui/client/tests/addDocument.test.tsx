import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockPostJson = vi.fn();
vi.mock('../src/io/entrypoint', () => ({
  postJson: (...args: any[]) => mockPostJson(...args),
  getJson: vi.fn(),
}));

const { AddDocumentButton } = await import('../src/components/project/AddDocumentButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');
const { useToolkitStore } = await import('../src/stores/useToolkitStore');

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
    currProjectName: 'TestProject',
    currProject: { name: 'TestProject', documents: [configDoc] },
    editedDocs: {},
  });
  useProjectListStore.setState({ projectNames: [{ name: 'TestProject' }] });
  useToolkitStore.setState({
    toolkits: [
      { toolkit: 'LSM', cls: 'lsm.cls' },
      { toolkit: 'GIS', cls: 'gis.cls' },
    ],
  });
});

const answerWith = (documents: any[], addedDocId = 'new1') => {
  mockPostJson.mockResolvedValueOnce({
    document: documents[documents.length - 1],
    project: { name: 'TestProject', documents },
    addedDocId,
  });
};

const madeDoc = (oid: string, desc: any, extra: any = {}) => ({
  _id: { $oid: oid },
  desc,
  type: 'T',
  _cls: 'Metadata.Measurements',
  resource: '',
  dataFormat: 'string',
  projectName: 'TestProject',
  ...extra,
});

const openDialog = async () => {
  render(<AddDocumentButton />);
  await act(async () => {
    fireEvent.click(screen.getByRole('button', { name: /add document/i }));
  });
  return await screen.findByRole('dialog');
};

const submit = async (dialog: HTMLElement) => {
  await act(async () => {
    fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
  });
};

const sentBody = () => {
  return mockPostJson.mock.calls[0][1];
};

describe('AddDocumentButton', () => {
  it('opens dialog when clicked', async () => {
    expect(await openDialog()).toBeTruthy();
  });

  it('calls document.create with the name, resource and collection', async () => {
    answerWith([madeDoc('new1', { datasourceName: 'MyDoc' }, { resource: 'myres' })]);
    const dialog = await openDialog();

    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'MyDoc' } });
    fireEvent.change(within(dialog).getByRole('textbox', { name: /resource/i }), { target: { value: 'myres' } });
    await submit(dialog);

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(mockPostJson.mock.calls[0][0]).toBe('/document/create');
    expect(sentBody()).toEqual({
      projectName: 'TestProject',
      collection: 'Measurements',
      resource: 'myres',
      desc: { datasourceName: 'MyDoc' },
    });
  });

  it('creates an agent document with effects resource', async () => {
    answerWith([madeDoc('ag1', { datasourceName: 'Agent1' })]);
    const dialog = await openDialog();

    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'Agent1' } });
    fireEvent.click(within(dialog).getByRole('button', { name: 'Agent' }));
    await submit(dialog);

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(sentBody().resource).toEqual({ effects: {} });
    expect(sentBody().type).toBe('ToolkitDataSource');
    expect(sentBody().dataFormat).toBe('JSON_dict');
  });

  it('creates a workflow document with an empty workflow block', async () => {
    answerWith([madeDoc('wf1', { datasourceName: 'Flow1' })]);
    const dialog = await openDialog();

    fireEvent.click(within(dialog).getByRole('button', { name: 'Workflow' }));
    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'Flow1' } });
    await submit(dialog);

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(sentBody().collection).toBe('Simulations');
    expect(sentBody().type).toBe('hermesWorkflow');
    expect(sentBody().resource).toBe('/tmp/testproject/Flow1.json');
    expect(sentBody().desc.workflowName).toBe('Flow1');
    expect(sentBody().desc.workflow).toEqual({ workflow: { solver: '', nodeList: [], nodes: {} } });
    expect(sentBody().desc.parameters).toEqual({});
    // No file is written for a workflow.
    expect(sentBody().notebook).toBeUndefined();
  });

  it('updates the store after successful creation', async () => {
    answerWith([madeDoc('n1', { datasourceName: 'NewDoc' })]);
    const dialog = await openDialog();

    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'NewDoc' } });
    await submit(dialog);

    await waitFor(() => {
      const project = useProjectStore.getState().currProject;
      expect(project?.documents).toHaveLength(1);
      expect(project?.documents[0].desc.datasourceName).toBe('NewDoc');
    });
  });

  it('hands the new document id to the caller', async () => {
    answerWith([madeDoc('n1', { datasourceName: 'NewDoc' })], 'n1');
    const onDocumentCreated = vi.fn();
    render(<AddDocumentButton onDocumentCreated={onDocumentCreated} />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /add document/i }));
    });
    const dialog = await screen.findByRole('dialog');

    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'NewDoc' } });
    await submit(dialog);

    await waitFor(() => {
      expect(onDocumentCreated).toHaveBeenCalledWith('n1');
    });
  });

  it('does not update store on error', async () => {
    mockPostJson.mockResolvedValueOnce(null);
    const dialog = await openDialog();

    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'Fail' } });
    await submit(dialog);

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(useProjectStore.getState().currProject?.documents).toEqual([configDoc]);
  });

  it('handles a delayed answer', async () => {
    let resolveCall!: (v: any) => void;
    mockPostJson.mockImplementationOnce(() => new Promise(r => { resolveCall = r; }));
    const dialog = await openDialog();

    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'Delayed' } });
    await submit(dialog);

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(useProjectStore.getState().currProject?.documents).toEqual([configDoc]);

    const documents = [madeDoc('d1', { datasourceName: 'Delayed' })];
    await act(async () => {
      resolveCall({ document: documents[0], project: { name: 'TestProject', documents }, addedDocId: 'd1' });
    });

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.documents).toHaveLength(1);
    });
  });

  it('creates a notebook document and asks the server for the file', async () => {
    answerWith([configDoc, madeDoc('nb1', { datasourceName: 'MyNotebook' })], 'nb1');
    const dialog = await openDialog();

    fireEvent.click(within(dialog).getByRole('button', { name: 'Notebook' }));
    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'MyNotebook' } });
    await submit(dialog);

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(sentBody().collection).toBe('Cache');
    expect(sentBody().type).toBe('notebook');
    expect(sentBody().resource).toBe('/tmp/testproject/notebooks/MyNotebook.ipynb');
    expect(sentBody().notebook).toEqual({ toolkitNames: [] });
  });

  it('shows resource field with auto-generated path when notebook is selected', async () => {
    const dialog = await openDialog();
    fireEvent.click(within(dialog).getByRole('button', { name: 'Notebook' }));

    const nameInput = within(dialog).getByRole('textbox', { name: /^name$/i });
    fireEvent.change(nameInput, { target: { value: 'TestNB' } });

    const resourceInput = within(dialog).getByRole('textbox', { name: /resource/i });
    expect(resourceInput).toHaveProperty('disabled', true);
    expect((resourceInput as HTMLInputElement).value).toBe('/tmp/testproject/notebooks/TestNB.ipynb');
  });

  it('shows helper text when notebook is selected', async () => {
    const dialog = await openDialog();
    fireEvent.click(within(dialog).getByRole('button', { name: 'Notebook' }));

    expect(within(dialog).getByText(/already exists at this path/i)).toBeTruthy();
  });

  it('hides class/toolkit when notebook kind is selected', async () => {
    const dialog = await openDialog();
    expect(within(dialog).queryByRole('combobox', { name: /toolkit/i })).toBeTruthy();

    fireEvent.click(within(dialog).getByRole('button', { name: 'Notebook' }));

    expect(within(dialog).queryByRole('combobox', { name: /toolkit/i })).toBeNull();
  });

  it('includes toolkit in desc when selected', async () => {
    answerWith([madeDoc('n1', { datasourceName: 'WithTK', toolkit: 'LSM' })]);
    render(<AddDocumentButton toolkit={{ toolkit: 'LSM', cls: 'lsm.cls' }} />);
    await act(async () => {
      fireEvent.click(screen.getByRole('button', { name: /add document/i }));
    });
    const dialog = await screen.findByRole('dialog');

    fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), { target: { value: 'WithTK' } });
    await submit(dialog);

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(sentBody().desc.toolkit).toBe('LSM');
  });

});
