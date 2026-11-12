import { describe, it, expect, vi, beforeEach } from 'vitest';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

import { mockApiCall } from './helpers/mockApiCall';

let mockPostJson: ReturnType<typeof vi.fn>;

const mockPushInfo = vi.fn();
vi.mock('../src/io/snackbar', () => ({
  pushInfo: (...args: any[]) => mockPushInfo(...args),
}));

const { detectNotebooks } = await import('../src/io/detectNotebooks');
const { useProjectStore } = await import('../src/stores/useProjectStore');

const configDoc = {
  _id: { $oid: 'cfg1' },
  desc: { datasourceName: 'config', filesDirectory: '/tmp/testproject' },
  type: 'TestProject__config__',
  _cls: 'Cache',
  resource: '',
  dataFormat: 'string',
  projectName: 'TestProject',
};

const newNotebook = {
  _id: { $oid: 'nb1' },
  desc: { datasourceName: 'Found' },
  type: 'notebook',
  _cls: 'Cache',
  resource: '/tmp/testproject/Found.ipynb',
  dataFormat: 'JSON_dict',
  projectName: 'TestProject',
};

const project = { name: 'TestProject', documents: [configDoc] };

beforeEach(() => {
  vi.clearAllMocks();
  mockPostJson = mockApiCall().postJson;
  useProjectStore.setState({ currProjectName: 'TestProject', currProject: project });
});

describe('detectNotebooks', () => {
  it('sends the project and its files directory to project.detectNotebooks', async () => {
    mockPostJson.mockResolvedValueOnce({ project, addedNotebooks: 0 });

    await detectNotebooks(project);

    expect(mockPostJson).toHaveBeenCalledTimes(1);
    const [path, body] = mockPostJson.mock.calls[0];
    expect(path).toBe('/project/detectNotebooks');
    expect(body).toEqual({ projectName: 'TestProject', filesDirectory: '/tmp/testproject' });
  });

  it('updates the store and reports the count when notebooks are found', async () => {
    mockPostJson.mockResolvedValueOnce({
      project: { name: 'TestProject', documents: [configDoc, newNotebook] },
      addedNotebooks: 1,
    });

    await detectNotebooks(project);

    expect(useProjectStore.getState().currProject?.documents).toHaveLength(2);
    expect(mockPushInfo).toHaveBeenCalledWith('Detected 1 new notebook');
  });

  it('says nothing when no new notebooks are found', async () => {
    mockPostJson.mockResolvedValueOnce({ project, addedNotebooks: 0 });

    await detectNotebooks(project);

    expect(mockPushInfo).not.toHaveBeenCalled();
  });

  it('does not touch the store when the user switched to another project', async () => {
    useProjectStore.setState({ currProjectName: 'Other', currProject: null });
    mockPostJson.mockResolvedValueOnce({
      project: { name: 'TestProject', documents: [configDoc, newNotebook] },
      addedNotebooks: 1,
    });

    await detectNotebooks(project);

    expect(useProjectStore.getState().currProject).toBeNull();
    expect(mockPushInfo).not.toHaveBeenCalled();
  });

  it('does not update the store when the call failed', async () => {
    mockPostJson.mockResolvedValueOnce(null);

    await detectNotebooks(project);

    expect(useProjectStore.getState().currProject?.documents).toEqual([configDoc]);
    expect(mockPushInfo).not.toHaveBeenCalled();
  });

  it('does not call the server when the project has no files directory', async () => {
    await detectNotebooks({ name: 'TestProject', documents: [] });

    expect(mockPostJson).not.toHaveBeenCalled();
  });
});
