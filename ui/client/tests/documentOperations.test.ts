import { describe, it, expect, vi, beforeEach } from 'vitest';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockPostJson = vi.fn();
vi.mock('../src/io/entrypoint', () => ({
  postJson: (...args: any[]) => mockPostJson(...args),
  getJson: vi.fn(),
}));

const { updateDocument } = await import('../src/io/FetchDocument');
const { deleteDocuments } = await import('../src/io/deleteDocuments');
const { duplicateDocuments, freeCopyName } = await import('../src/io/duplicateDocuments');
const { useProjectStore } = await import('../src/stores/useProjectStore');

const doc = (oid: string, name: string, extra: any = {}) => ({
  _id: { $oid: oid },
  _cls: 'Metadata.Measurements',
  projectName: 'P',
  desc: { datasourceName: name },
  type: 'T',
  resource: `res-${oid}`,
  dataFormat: 'string',
  ...extra,
});

const setProject = (documents: any[]) => {
  useProjectStore.setState({
    currProjectName: 'P',
    currProject: { name: 'P', documents },
    editedDocs: {},
  });
};

beforeEach(() => {
  vi.clearAllMocks();
  setProject([]);
});

const sent = (call = 0) => {
  return mockPostJson.mock.calls[call][1];
};

describe('updateDocument', () => {
  const stored = doc('doc1', 'old');

  it('sends the document id and the whole document as fields', async () => {
    mockPostJson.mockResolvedValueOnce({
      document: { ...stored, resource: 'res2' },
      project: { name: 'P', documents: [] },
    });

    await updateDocument({ ...stored, resource: 'res2' } as any);

    expect(mockPostJson.mock.calls[0][0]).toBe('/document/update');
    expect(sent().docId).toBe('doc1');
    expect(sent().fields.resource).toBe('res2');
    expect(sent().fields.type).toBe('T');
    expect(sent().fields.desc).toEqual({ datasourceName: 'old' });
  });

  it('never sends _id, _cls or projectName', async () => {
    mockPostJson.mockResolvedValueOnce({ document: stored, project: { name: 'P', documents: [] } });

    await updateDocument({ ...stored, projectName: 'Changed' } as any);

    expect(sent().fields._id).toBeUndefined();
    expect(sent().fields._cls).toBeUndefined();
    expect(sent().fields.projectName).toBeUndefined();
  });

  it('carries null and booleans as they are', async () => {
    mockPostJson.mockResolvedValueOnce({ document: stored, project: { name: 'P', documents: [] } });

    await updateDocument({
      ...stored,
      resource: null,
      desc: { datasourceName: 'x', flag: true, disabled: false },
    } as any);

    expect(sent().fields.resource).toBeNull();
    expect(sent().fields.desc.flag).toBe(true);
    expect(sent().fields.desc.disabled).toBe(false);
  });

  it('answers the document and the project', async () => {
    const saved = { ...stored, resource: 'res2' };
    mockPostJson.mockResolvedValueOnce({ document: saved, project: { name: 'P', documents: [saved] } });

    const answer = await updateDocument(saved as any);

    expect(answer?.document).toEqual(saved);
    expect(answer?.project.documents).toEqual([saved]);
  });

  it('answers null when the call fails', async () => {
    mockPostJson.mockResolvedValueOnce(null);

    expect(await updateDocument(stored as any)).toBeNull();
  });
});

describe('deleteDocuments', () => {
  it('sends every id in one call and stores the project it answers', async () => {
    setProject([doc('doc1', 'a'), doc('doc2', 'b')]);
    mockPostJson.mockResolvedValueOnce({ project: { name: 'P', documents: [] } });

    expect(await deleteDocuments(['doc1', 'doc2'])).toBe(true);

    expect(mockPostJson).toHaveBeenCalledTimes(1);
    expect(mockPostJson.mock.calls[0][0]).toBe('/document/delete');
    expect(sent()).toEqual({ docIds: ['doc1', 'doc2'] });
    expect(useProjectStore.getState().currProject?.documents).toEqual([]);
  });

  it('says so and keeps the project when the call fails', async () => {
    const documents = [doc('doc1', 'a')];
    setProject(documents);
    mockPostJson.mockResolvedValueOnce(null);

    expect(await deleteDocuments(['doc1'])).toBe(false);
    expect(useProjectStore.getState().currProject?.documents).toEqual(documents);
  });
});

describe('freeCopyName', () => {
  it('is "<name> copy" when nothing takes it', () => {
    expect(freeCopyName('MyDoc', ['MyDoc'])).toBe('MyDoc copy');
  });

  it('counts up past the copies already there', () => {
    expect(freeCopyName('MyDoc', ['MyDoc', 'MyDoc copy'])).toBe('MyDoc copy 2');
    expect(freeCopyName('MyDoc', ['MyDoc', 'MyDoc copy', 'MyDoc copy 2'])).toBe('MyDoc copy 3');
  });

  it('skips a gap it does not need to fill', () => {
    expect(freeCopyName('MyDoc', ['MyDoc copy', 'MyDoc copy 3'])).toBe('MyDoc copy 2');
  });
});

describe('duplicateDocuments', () => {
  const answerWith = (documents: any[]) => {
    mockPostJson.mockResolvedValueOnce({
      document: documents[documents.length - 1],
      project: { name: 'P', documents },
      addedDocId: 'copy1',
    });
  };

  it('copies the stored fields into the document collection', async () => {
    const original = doc('doc1', 'MyDoc', { type: 'ToolkitDataSource', dataFormat: 'JSON_dict' });
    setProject([original]);
    answerWith([original, doc('copy1', 'MyDoc copy')]);

    expect(await duplicateDocuments(['doc1'])).toBe(true);

    expect(mockPostJson.mock.calls[0][0]).toBe('/document/create');
    expect(sent()).toEqual({
      projectName: 'P',
      collection: 'Measurements',
      type: 'ToolkitDataSource',
      resource: 'res-doc1',
      dataFormat: 'JSON_dict',
      desc: { datasourceName: 'MyDoc copy' },
    });
    expect(useProjectStore.getState().currProject?.documents).toHaveLength(2);
  });

  it('reads nothing back: one call per document', async () => {
    setProject([doc('doc1', 'A'), doc('doc2', 'B')]);
    answerWith([doc('copy1', 'A copy')]);
    answerWith([doc('copy2', 'B copy')]);

    await duplicateDocuments(['doc1', 'doc2']);

    expect(mockPostJson).toHaveBeenCalledTimes(2);
    expect(sent(0).desc.datasourceName).toBe('A copy');
    expect(sent(1).desc.datasourceName).toBe('B copy');
  });

  it('gives two copies of one document different names', async () => {
    setProject([doc('doc1', 'A')]);
    answerWith([doc('copy1', 'A copy')]);
    answerWith([doc('copy2', 'A copy 2')]);

    await duplicateDocuments(['doc1', 'doc1']);

    expect(sent(0).desc.datasourceName).toBe('A copy');
    expect(sent(1).desc.datasourceName).toBe('A copy 2');
  });

  it('maps the document class to its collection', async () => {
    setProject([doc('doc1', 'A', { _cls: 'Metadata.Simulations' })]);
    answerWith([doc('copy1', 'A copy')]);

    await duplicateDocuments(['doc1']);

    expect(sent().collection).toBe('Simulations');
  });

  it('stops on a class it cannot place', async () => {
    setProject([doc('doc1', 'A', { _cls: 'Metadata' })]);

    expect(await duplicateDocuments(['doc1'])).toBe(false);
    expect(mockPostJson).not.toHaveBeenCalled();
  });

  it('stops at the first failed call', async () => {
    setProject([doc('doc1', 'A'), doc('doc2', 'B')]);
    mockPostJson.mockResolvedValueOnce(null);

    expect(await duplicateDocuments(['doc1', 'doc2'])).toBe(false);
    expect(mockPostJson).toHaveBeenCalledTimes(1);
  });
});
