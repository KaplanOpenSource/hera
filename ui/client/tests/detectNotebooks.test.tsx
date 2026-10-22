import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, fireEvent, waitFor, cleanup, act } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

import { mockApiCall } from './helpers/mockApiCall';

let mockPostJson: ReturnType<typeof vi.fn>;

const mockPushInfo = vi.fn();
vi.mock('../src/io/snackbar', () => ({
  pushInfo: (...args: any[]) => mockPushInfo(...args),
}));

const { DetectNotebooksButton } = await import('../src/components/project/DetectNotebooksButton');
const { useProjectStore } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');

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
  mockPostJson = mockApiCall().postJson;
  useProjectStore.setState({
    currProjectName: 'TestProject',
    currProject: { name: 'TestProject', documents: [configDoc] },
  });
  useProjectListStore.setState({ projectNames: [{ name: 'TestProject' }] });
});

const clickDetect = async () => {
  await act(async () => {
    fireEvent.click(screen.getByRole('button', { name: /detect notebooks/i }));
  });
};

describe('DetectNotebooksButton', () => {
  it('sends the project and its files directory to project.detectNotebooks', async () => {
    mockPostJson.mockResolvedValueOnce({
      project: { name: 'TestProject', documents: [configDoc] },
      addedNotebooks: 0,
    });

    render(<DetectNotebooksButton />);
    await clickDetect();

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    const [path, body] = mockPostJson.mock.calls[0];
    expect(path).toBe('/project/detectNotebooks');
    expect(body).toEqual({ projectName: 'TestProject', filesDirectory: '/tmp/testproject' });
  });

  it('updates the store and reports the count when notebooks are found', async () => {
    const newNotebook = { _id: { $oid: 'nb1' }, desc: { datasourceName: 'Found' }, type: 'notebook', _cls: 'Cache', resource: '/tmp/testproject/Found.ipynb', dataFormat: 'JSON_dict', projectName: 'TestProject' };
    mockPostJson.mockResolvedValueOnce({
      project: { name: 'TestProject', documents: [configDoc, newNotebook] },
      addedNotebooks: 1,
    });

    const onDetected = vi.fn();
    render(<DetectNotebooksButton onDetected={onDetected} />);
    await clickDetect();

    await waitFor(() => {
      expect(useProjectStore.getState().currProject?.documents).toHaveLength(2);
      expect(mockPushInfo).toHaveBeenCalledWith('Detected 1 new notebook');
      expect(onDetected).toHaveBeenCalled();
    });
  });

  it('reports when no new notebooks are found', async () => {
    mockPostJson.mockResolvedValueOnce({
      project: { name: 'TestProject', documents: [configDoc] },
      addedNotebooks: 0,
    });

    render(<DetectNotebooksButton />);
    await clickDetect();

    await waitFor(() => {
      expect(mockPushInfo).toHaveBeenCalledWith('No new notebooks found');
    });
  });

  it('does not update the store when the call failed', async () => {
    mockPostJson.mockResolvedValueOnce(null);

    render(<DetectNotebooksButton />);
    await clickDetect();

    await waitFor(() => {
      expect(mockPostJson).toHaveBeenCalledTimes(1);
    });
    expect(useProjectStore.getState().currProject?.documents).toEqual([configDoc]);
    expect(mockPushInfo).not.toHaveBeenCalled();
  });
});
