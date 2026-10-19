import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { act, cleanup, render } from '@testing-library/react';
import { mockApiCall } from './helpers/mockApiCall';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const { useShownDoc } = await import('../src/components/details/useShownDoc');
const { ProjectObj } = await import('../src/objects/ProjectObj');
const { useProjectStore } = await import('../src/stores/useProjectStore');

afterEach(() => { cleanup(); });

let mockGetJson: ReturnType<typeof vi.fn>;
let mockPostJson: ReturnType<typeof vi.fn>;

const doc = (resource: string) => ({
  _id: { $oid: 'doc1' },
  _cls: 'Metadata.Measurements',
  projectName: 'P',
  desc: { datasourceName: 'MyDoc' },
  type: 'T',
  resource,
  dataFormat: 'string',
});

const project = (resource: string) => ({ name: 'P', documents: [doc(resource)] });

beforeEach(() => {
  vi.clearAllMocks();
  const mocked = mockApiCall();
  mockGetJson = mocked.getJson;
  mockPostJson = mocked.postJson;
  useProjectStore.setState({
    currProjectName: 'P',
    currProject: project('old') as any,
    editedDocs: {},
  });
});

// Drives the hook from a throwaway component.
const mountHook = () => {
  const held: { save?: () => Promise<void>, edit?: (resource: string) => void } = {};
  const Probe = () => {
    const docObj = new ProjectObj(project('old') as any).allDocuments[0];
    const { setShownDoc, saveShownDoc } = useShownDoc(docObj);
    held.save = saveShownDoc;
    held.edit = (resource: string) => {
      return setShownDoc({ ...docObj.data, resource } as any);
    };
    return null;
  };
  render(<Probe />);
  return held;
};

describe('saving a document', () => {
  it('is one call: the answer carries the project', async () => {
    const held = mountHook();
    await act(async () => { held.edit!('new'); });

    mockPostJson.mockResolvedValueOnce({
      document: doc('new'),
      project: project('new'),
    });
    await act(async () => { await held.save!(); });

    expect(mockPostJson).toHaveBeenCalledTimes(1);
    expect(mockPostJson.mock.calls[0][0]).toBe('/document/update');
    // No reload after it.
    expect(mockGetJson).not.toHaveBeenCalled();
    expect(useProjectStore.getState().currProject?.documents[0].resource).toBe('new');
    expect(useProjectStore.getState().editedDocs['doc1']).toBeUndefined();
  });

  it('keeps the edits when the save fails', async () => {
    const held = mountHook();
    await act(async () => { held.edit!('new'); });

    mockPostJson.mockResolvedValueOnce(null);
    await act(async () => { await held.save!(); });

    expect(useProjectStore.getState().editedDocs['doc1']?.resource).toBe('new');
    expect(useProjectStore.getState().currProject?.documents[0].resource).toBe('old');
  });

  it('calls nothing when there is nothing to save', async () => {
    const held = mountHook();

    await act(async () => { await held.save!(); });

    expect(mockPostJson).not.toHaveBeenCalled();
  });
});
