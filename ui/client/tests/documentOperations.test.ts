import { describe, it, expect, vi, beforeEach } from 'vitest';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetch = vi.fn();
global.fetch = mockFetch;

const { fetchDocument, updateDocument } = await import('../src/io/FetchDocument');

beforeEach(() => {
  vi.clearAllMocks();
});

describe('fetchDocument', () => {
  it('returns document data on success', async () => {
    const doc = { _id: { $oid: 'abc' }, type: 'T', desc: {} };
    mockFetch.mockResolvedValueOnce({ json: () => Promise.resolve(doc) });

    const result = await fetchDocument('abc');

    expect(result).toEqual(doc);
    const body = JSON.parse(mockFetch.mock.calls[0][1].body);
    expect(body.code).toContain("All.getDocumentByID('abc')");
  });

  it('returns undefined when data is falsy', async () => {
    mockFetch.mockResolvedValueOnce({ json: () => Promise.resolve(null) });
    const result = await fetchDocument('abc');
    expect(result).toBeUndefined();
  });

  it('works with delayed response', async () => {
    const doc = { _id: { $oid: 'abc' }, type: 'T', desc: {} };
    mockFetch.mockImplementationOnce(() =>
      new Promise(resolve => setTimeout(() =>
        resolve({ json: () => Promise.resolve(doc) }), 50
      ))
    );

    const result = await fetchDocument('abc');

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
    const updated = { ...newDoc };
    mockFetch.mockResolvedValueOnce({ json: () => Promise.resolve(updated) });

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
    mockFetch.mockResolvedValueOnce({ json: () => Promise.resolve(newDoc) });

    await updateDocument(newDoc, prevDoc);

    const body = JSON.parse(mockFetch.mock.calls[0][1].body);
    expect(body.code).not.toContain('doc._id');
    expect(body.code).not.toContain('doc._cls');
    expect(body.code).not.toContain('doc.projectName');
    expect(body.code).toContain('doc.resource');
  });

  it('converts null to Python None', async () => {
    const newDoc = { ...prevDoc, resource: null };
    mockFetch.mockResolvedValueOnce({ json: () => Promise.resolve(newDoc) });

    await updateDocument(newDoc, prevDoc);

    const body = JSON.parse(mockFetch.mock.calls[0][1].body);
    expect(body.code).toContain('doc.resource = None');
  });

  it('handles nested object changes', async () => {
    const newDoc = { ...prevDoc, desc: { datasourceName: 'new' } };
    mockFetch.mockResolvedValueOnce({ json: () => Promise.resolve(newDoc) });

    await updateDocument(newDoc, prevDoc);

    const body = JSON.parse(mockFetch.mock.calls[0][1].body);
    expect(body.code).toContain('doc.desc');
    expect(body.code).toContain('"datasourceName"');
  });

  it('calls save and refetches', async () => {
    const newDoc = { ...prevDoc, resource: 'changed' };
    mockFetch.mockResolvedValueOnce({ json: () => Promise.resolve(newDoc) });

    await updateDocument(newDoc, prevDoc);

    const body = JSON.parse(mockFetch.mock.calls[0][1].body);
    expect(body.code).toContain('doc.save()');
    expect(body.code).toContain("All.getDocumentByID('doc1')");
  });

  it('works with delayed response', async () => {
    const newDoc = { ...prevDoc, resource: 'delayed' };
    mockFetch.mockImplementationOnce(() =>
      new Promise(resolve => setTimeout(() =>
        resolve({ json: () => Promise.resolve(newDoc) }), 50
      ))
    );

    const result = await updateDocument(newDoc, prevDoc);

    expect(result).toEqual(newDoc);
  });
});
