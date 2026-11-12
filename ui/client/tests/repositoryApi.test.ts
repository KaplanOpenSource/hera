import { describe, it, expect, vi, beforeEach } from 'vitest';
import { mockApiCall } from './helpers/mockApiCall';

const { repositoryApi } = await import('../src/api/api');

let getJson: ReturnType<typeof vi.fn>;
let postJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson, postJson } = mockApiCall());
});

describe('repositoryApi.get', () => {
  it('sends the path and returns both trees', async () => {
    getJson.mockResolvedValueOnce({ json: { a: 1 }, resolved: { a: 2 }, extra: 3 });

    expect(await repositoryApi.get('/repos/r.json')).toEqual({ json: { a: 1 }, resolved: { a: 2 } });
    expect(getJson).toHaveBeenCalledWith('/repository/get', { path: '/repos/r.json' }, { label: 'repo JSON' });
  });

  it('returns null when the call failed', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await repositoryApi.get('/repos/r.json')).toBeNull();
  });
});

describe('repositoryApi.registered', () => {
  it('returns the repositories', async () => {
    getJson.mockResolvedValueOnce({ repositories: [{ name: 'r1' }] });

    expect(await repositoryApi.registered()).toEqual([{ name: 'r1' }]);
    expect(getJson).toHaveBeenCalledWith('/repository/registered', {}, { label: 'registered repositories' });
  });

  it('returns an empty list when the answer has none', async () => {
    getJson.mockResolvedValueOnce({});

    expect(await repositoryApi.registered()).toEqual([]);
  });

  it('returns an empty list when the call failed', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await repositoryApi.registered()).toEqual([]);
  });
});

describe('repositoryApi.jsonsInFolder', () => {
  it('asks for the file list only by default', async () => {
    getJson.mockResolvedValueOnce({ files: ['/f/a.json'] });

    expect(await repositoryApi.jsonsInFolder('/f')).toEqual({ files: ['/f/a.json'], repoJsons: {} });
    expect(getJson).toHaveBeenCalledWith('/repository/jsonsInFolder',
      { folder: '/f', withContents: 'false' },
      { label: 'central repo files' });
  });

  it('can ask for the contents too', async () => {
    getJson.mockResolvedValueOnce({ files: ['/f/a.json'], repoJsons: { '/f/a.json': { a: 1 } } });

    const answer = await repositoryApi.jsonsInFolder('/f', true);

    expect(answer).toEqual({ files: ['/f/a.json'], repoJsons: { '/f/a.json': { a: 1 } } });
    expect(getJson.mock.calls[0][1]).toEqual({ folder: '/f', withContents: 'true' });
  });

  it('fills in empty lists when the answer has none', async () => {
    getJson.mockResolvedValueOnce({});

    expect(await repositoryApi.jsonsInFolder('/f')).toEqual({ files: [], repoJsons: {} });
  });

  it('returns null when the call failed', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await repositoryApi.jsonsInFolder('/f')).toBeNull();
  });
});

describe('repositoryApi.loadIntoProject', () => {
  it('sends the body as it is, with the label', async () => {
    postJson.mockResolvedValueOnce({ project: { projectName: 'A' }, added: 1, problems: [] });
    const body = { projectName: 'A', path: '/repos/r.json', overwrite: false };

    expect(await repositoryApi.loadIntoProject(body, 'load r')).toEqual(
      { project: { projectName: 'A' }, added: 1, problems: [] });
    expect(postJson).toHaveBeenCalledWith('/repository/loadIntoProject', body, { label: 'load r' });
  });

  it('returns null when the call failed', async () => {
    postJson.mockResolvedValueOnce(null);

    expect(await repositoryApi.loadIntoProject({}, 'load')).toBeNull();
  });
});
