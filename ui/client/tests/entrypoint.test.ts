import { describe, it, expect, vi, beforeEach } from 'vitest';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockPushError = vi.fn();
const mockPushRunning = vi.fn().mockReturnValue('key1');
const mockDismiss = vi.fn();
vi.mock('../src/io/snackbar', () => ({
  pushError: (...args: any[]) => mockPushError(...args),
  pushRunning: (...args: any[]) => mockPushRunning(...args),
  dismiss: (...args: any[]) => mockDismiss(...args),
}));

const mockFetch = vi.fn();
global.fetch = mockFetch;

const { getJson, postJson } = await import('../src/io/entrypoint');

const jsonAnswer = (data: any, ok = true, status = 200) => ({
  ok,
  status,
  headers: { get: () => 'application/json' },
  json: () => Promise.resolve(data),
});

beforeEach(() => {
  vi.clearAllMocks();
});

describe('getJson', () => {
  it('returns the answer and shows a running toast while it waits', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({ projects: [] }));

    expect(await getJson('/project/names', {}, { label: 'names' })).toEqual({ projects: [] });
    expect(mockPushRunning).toHaveBeenCalledWith('names');
    expect(mockDismiss).toHaveBeenCalledWith('key1');
  });

  it('puts the parameters in the query, encoded', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({}));

    await getJson('/project/documents', { projectName: 'My Project' }, { label: 'docs' });

    expect(mockFetch.mock.calls[0][0]).toBe('http://test/project/documents?projectName=My+Project');
  });

  it('shows no toast at all when silent', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({}));

    await getJson('/project/names', {}, { label: 'names', silent: true });

    expect(mockPushRunning).not.toHaveBeenCalled();
    expect(mockDismiss).not.toHaveBeenCalled();
  });

  it('reports the server message and returns null on a rejected call', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({ detail: 'no such project: Nope' }, false, 404));

    expect(await getJson('/project/documents', { projectName: 'Nope' }, { label: 'docs' })).toBeNull();
    expect(mockPushError).toHaveBeenCalledWith('docs: no such project: Nope');
  });

  it('stays quiet while the server is warming up', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({ detail: 'WARMING_UP' }, false, 503));

    expect(await getJson('/project/names', {}, { label: 'names' })).toBeNull();
    expect(mockPushError).not.toHaveBeenCalled();
  });

  it('returns null when the server has no such route and answers HTML', async () => {
    mockFetch.mockResolvedValueOnce({
      ok: true,
      status: 200,
      headers: { get: () => 'text/html' },
      json: () => Promise.reject(new Error('not json')),
    });

    expect(await getJson('/project/names', {}, { label: 'names' })).toBeNull();
    expect(mockPushError).toHaveBeenCalled();
  });

  it('returns null when the server is not there', async () => {
    mockFetch.mockRejectedValueOnce(new Error('Failed to fetch'));

    expect(await getJson('/project/names', {}, { label: 'names' })).toBeNull();
    expect(mockPushError).toHaveBeenCalledWith('names: network error: Failed to fetch');
    expect(mockDismiss).toHaveBeenCalledWith('key1');
  });
});

describe('postJson', () => {
  it('sends the body as JSON', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({ project: null }));

    await postJson('/project/delete', { projectName: 'Alpha', deleteFiles: true }, { label: 'delete' });

    const [url, init] = mockFetch.mock.calls[0];
    expect(url).toBe('http://test/project/delete');
    expect(init.method).toBe('POST');
    expect(init.headers['Content-Type']).toBe('application/json');
    expect(JSON.parse(init.body)).toEqual({ projectName: 'Alpha', deleteFiles: true });
  });
});
