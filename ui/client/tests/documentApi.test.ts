import { describe, it, expect, vi, beforeEach } from 'vitest';
import { mockApiCall } from './helpers/mockApiCall';
import { Collection } from '../src/shared/types';

const { documentApi } = await import('../src/api/api');

let postJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  ({ postJson } = mockApiCall());
});

describe('documentApi.create', () => {
  it('sends the new document with the label', async () => {
    postJson.mockResolvedValueOnce({ project: { projectName: 'A' }, addedDocId: 'd1' });
    const newDocument = { projectName: 'A', collection: 'Cache' as Collection, type: 'T', desc: { a: 1 } };

    expect(await documentApi.create(newDocument, 'add T')).toEqual({ project: { projectName: 'A' }, addedDocId: 'd1' });
    expect(postJson).toHaveBeenCalledWith('/document/create', newDocument, { label: 'add T' });
  });

  it('keeps only the project and the new id', async () => {
    postJson.mockResolvedValueOnce({ project: { projectName: 'A' }, addedDocId: 'd1', extra: 1 });

    const answer = await documentApi.create({ projectName: 'A', collection: 'Cache' as Collection }, 'add');

    expect(answer).toEqual({ project: { projectName: 'A' }, addedDocId: 'd1' });
  });

  it('returns null when the call failed', async () => {
    postJson.mockResolvedValueOnce(null);

    expect(await documentApi.create({ projectName: 'A', collection: 'Cache' as Collection }, 'add')).toBeNull();
  });
});

describe('documentApi.update', () => {
  const doc = {
    _id: { $oid: 'd1' },
    _cls: 'Metadata.Cache',
    projectName: 'A',
    type: 'T',
    desc: { a: 1 },
    resource: '/r',
  } as any;

  it('sends the id and every field but the ones that name the document', async () => {
    postJson.mockResolvedValueOnce({ document: doc, project: { projectName: 'A' } });

    await documentApi.update(doc);

    expect(postJson).toHaveBeenCalledWith('/document/update',
      { docId: 'd1', fields: { type: 'T', desc: { a: 1 }, resource: '/r' } },
      { label: 'update document' });
  });

  it('returns the saved document and the project', async () => {
    postJson.mockResolvedValueOnce({ document: doc, project: { projectName: 'A' }, extra: 1 });

    expect(await documentApi.update(doc)).toEqual({ document: doc, project: { projectName: 'A' } });
  });

  it('returns null when the call failed', async () => {
    postJson.mockResolvedValueOnce(null);

    expect(await documentApi.update(doc)).toBeNull();
  });
});

describe('documentApi.remove', () => {
  it('sends the ids and counts them in the label', async () => {
    postJson.mockResolvedValueOnce({ project: { projectName: 'A' } });

    expect(await documentApi.remove(['d1', 'd2'])).toEqual({ project: { projectName: 'A' } });
    expect(postJson).toHaveBeenCalledWith('/document/delete', { docIds: ['d1', 'd2'] },
      { label: 'delete 2 documents' });
  });

  it('returns null when the call failed', async () => {
    postJson.mockResolvedValueOnce(null);

    expect(await documentApi.remove(['d1'])).toBeNull();
  });
});
