import { describe, it, expect, vi, beforeEach } from 'vitest';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetch = vi.fn();
global.fetch = mockFetch;

const { ApiCall } = await import('../src/api/ApiCall');

const started = vi.fn().mockReturnValue('key1');
const finished = vi.fn();
const failed = vi.fn();
const reporter = { started, finished, failed };

// The protected methods are reachable at run time; a test subclass opens them.
class TestApi extends ApiCall {
  get(path: string, params: any, options: any) {
    return this.getJson(path, params, options);
  }
  post(path: string, body: any, options: any) {
    return this.postJson(path, body, options);
  }
}

const jsonAnswer = (data: any, ok = true, status = 200) => ({
  ok,
  status,
  headers: { get: () => 'application/json' },
  json: () => Promise.resolve(data),
});

let api: TestApi;

beforeEach(() => {
  vi.clearAllMocks();
  started.mockReturnValue('key1');
  api = new TestApi(reporter);
});

describe('getJson', () => {
  it('returns the answer and reports a running call while it waits', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({ projects: [] }));

    expect(await api.get('/project/names', {}, { label: 'names' })).toEqual({ projects: [] });
    expect(started).toHaveBeenCalledWith('names');
    expect(finished).toHaveBeenCalledWith('key1');
  });

  it('puts the parameters in the query, encoded', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({}));

    await api.get('/project/documents', { projectName: 'My Project' }, { label: 'docs' });

    expect(mockFetch.mock.calls[0][0]).toBe('http://test/project/documents?projectName=My+Project');
  });

  it('reports nothing at all when silent', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({}));

    await api.get('/project/names', {}, { label: 'names', silent: true });

    expect(started).not.toHaveBeenCalled();
    expect(finished).not.toHaveBeenCalled();
  });

  it('reports the server message and returns null on a rejected call', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({ detail: 'no such project: Nope' }, false, 404));

    expect(await api.get('/project/documents', { projectName: 'Nope' }, { label: 'docs' })).toBeNull();
    expect(failed).toHaveBeenCalledWith('docs: no such project: Nope');
  });

  it('stays quiet while the server is warming up', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({ detail: 'WARMING_UP' }, false, 503));

    expect(await api.get('/project/names', {}, { label: 'names' })).toBeNull();
    expect(failed).not.toHaveBeenCalled();
  });

  it('returns null when the server has no such route and answers HTML', async () => {
    mockFetch.mockResolvedValueOnce({
      ok: true,
      status: 200,
      headers: { get: () => 'text/html' },
      json: () => Promise.reject(new Error('not json')),
    });

    expect(await api.get('/project/names', {}, { label: 'names' })).toBeNull();
    expect(failed).toHaveBeenCalled();
  });

  it('returns null when the server is not there', async () => {
    mockFetch.mockRejectedValueOnce(new Error('Failed to fetch'));

    expect(await api.get('/project/names', {}, { label: 'names' })).toBeNull();
    expect(failed).toHaveBeenCalledWith('names: network error: Failed to fetch');
    expect(finished).toHaveBeenCalledWith('key1');
  });

  it('works with no reporter installed', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({ projects: [] }));
    const quiet = new TestApi();

    expect(await quiet.get('/project/names', {}, { label: 'names' })).toEqual({ projects: [] });
    expect(started).not.toHaveBeenCalled();
  });
});

describe('postJson', () => {
  it('sends the body as JSON', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({ project: null }));

    await api.post('/project/delete', { projectName: 'Alpha', deleteFiles: true }, { label: 'delete' });

    const [url, init] = mockFetch.mock.calls[0];
    expect(url).toBe('http://test/project/delete');
    expect(init.method).toBe('POST');
    expect(init.headers['Content-Type']).toBe('application/json');
    expect(JSON.parse(init.body)).toEqual({ projectName: 'Alpha', deleteFiles: true });
  });
});

describe('setReporter', () => {
  it('installs a reporter after the object was built', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({}));
    const late = new TestApi();
    late.setReporter(reporter);

    await late.get('/project/names', {}, { label: 'names' });

    expect(started).toHaveBeenCalledWith('names');
  });
});
