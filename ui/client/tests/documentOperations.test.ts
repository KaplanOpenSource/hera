import { describe, it, expect, vi, beforeEach } from 'vitest';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));
vi.mock('../src/io/snackbar', () => ({
  pushRunning: () => 'mock-key',
  pushError: () => 'mock-key',
  dismiss: () => {},
}));

const mockFetch = vi.fn();
global.fetch = mockFetch;

const { fetchDocument, updateDocument } = await import('../src/io/FetchDocument');

const mockResponse = (data: any) => ({
  ok: true,
  text: () => Promise.resolve(JSON.stringify({ data, problem: null })),
});

beforeEach(() => {
  vi.clearAllMocks();
});

describe('fetchDocument', () => {
  it('returns document data on success', async () => {
    const doc = { _id: { $oid: 'abc' }, type: 'T', desc: {} };
    mockFetch.mockResolvedValueOnce(mockResponse({ docData: doc }));

    const result = await fetchDocument('abc');

    expect(result).toEqual(doc);
    const body = JSON.parse(mockFetch.mock.calls[0][1].body);
    expect(body.code).toContain("All.getDocumentByID('abc')");
  });

  it('returns undefined when data is falsy', async () => {
    mockFetch.mockResolvedValueOnce(mockResponse(null));
    const result = await fetchDocument('abc');
    expect(result).toBeUndefined();
  });

  it('result not available while response is in flight', async () => {
    const doc = { _id: { $oid: 'abc' }, type: 'T', desc: {} };
    let resolve!: (v: any) => void;
    mockFetch.mockImplementationOnce(() =>
      new Promise(r => { resolve = r; })
    );

    let result: any;
    const promise = fetchDocument('abc').then(r => { result = r; });
    expect(result).toBeUndefined();

    resolve(mockResponse({ docData: doc }));
    await promise;

    expect(result).toEqual(doc);
  });
});

describe('updateDocument', () => {
  const prevDoc = {
    _id: { $oid: 'doc1' },
    _cls: 'Metadata.Cache',
    projectName: 'P',
    desc: { datasourceName: 'old' },
    type: 'T',
    resource: 'res1',
    dataFormat: 'string',
  };

  it('sends only changed fields', async () => {
    const newDoc = { ...prevDoc, resource: 'res2' };
    mockFetch.mockResolvedValueOnce(mockResponse({ docData: newDoc }));

    await updateDocument(newDoc, prevDoc);

    const body = JSON.parse(mockFetch.mock.calls[0][1].body);
    expect(body.code).toContain("doc.resource = \"res2\"");
    expect(body.code).not.toContain('doc.type');
    expect(body.code).not.toContain('doc.desc');
  });

  it('skips forbidden fields (_id, _cls, projectName)', async () => {
    const newDoc = {
      ...prevDoc,
      _id: { $oid: 'changed' },
      _cls: 'Changed',
      projectName: 'Changed',
      resource: 'new',
    };
    mockFetch.mockResolvedValueOnce(mockResponse({ docData: newDoc }));

    await updateDocument(newDoc, prevDoc);

    const body = JSON.parse(mockFetch.mock.calls[0][1].body);
    expect(body.code).not.toContain('doc._id');
    expect(body.code).not.toContain('doc._cls');
    expect(body.code).not.toContain('doc.projectName');
    expect(body.code).toContain('doc.resource');
  });

  it('converts null to Python None', async () => {
    const newDoc = { ...prevDoc, resource: null };
    mockFetch.mockResolvedValueOnce(mockResponse({ docData: newDoc }));

    await updateDocument(newDoc, prevDoc);

    const body = JSON.parse(mockFetch.mock.calls[0][1].body);
    expect(body.code).toContain('doc.resource = None');
  });

  it('converts boolean true/false to Python True/False', async () => {
    const newDoc = { ...prevDoc, desc: { datasourceName: 'changed', flag: true, disabled: false } };
    mockFetch.mockResolvedValueOnce(mockResponse({ docData: newDoc }));

    await updateDocument(newDoc, prevDoc);

    const body = JSON.parse(mockFetch.mock.calls[0][1].body);
    const code: string = body.code;
    expect(code).not.toMatch(/:true[,}]/);
    expect(code).not.toMatch(/:false[,}]/);
    expect(code).toMatch(/:True[,}]/);
    expect(code).toMatch(/:False[,}]/);
  });

  it('handles nested object changes', async () => {
    const newDoc = { ...prevDoc, desc: { datasourceName: 'new' } };
    mockFetch.mockResolvedValueOnce(mockResponse({ docData: newDoc }));

    await updateDocument(newDoc, prevDoc);

    const body = JSON.parse(mockFetch.mock.calls[0][1].body);
    expect(body.code).toContain('doc.desc');
    expect(body.code).toContain('"datasourceName"');
  });

  it('calls save and refetches', async () => {
    const newDoc = { ...prevDoc, resource: 'changed' };
    mockFetch.mockResolvedValueOnce(mockResponse({ docData: newDoc }));

    await updateDocument(newDoc, prevDoc);

    const body = JSON.parse(mockFetch.mock.calls[0][1].body);
    expect(body.code).toContain('doc.save()');
    expect(body.code).toContain("All.getDocumentByID('doc1')");
  });

  it('result not available while response is in flight', async () => {
    const newDoc = { ...prevDoc, resource: 'delayed' };
    let resolve!: (v: any) => void;
    mockFetch.mockImplementationOnce(() =>
      new Promise(r => { resolve = r; })
    );

    let result: any;
    const promise = updateDocument(newDoc, prevDoc).then(r => { result = r; });
    expect(result).toBeUndefined();

    resolve(mockResponse({ docData: newDoc }));
    await promise;

    expect(result).toEqual(newDoc);
  });
});
