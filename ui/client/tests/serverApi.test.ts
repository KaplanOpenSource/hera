import { describe, it, expect, vi, beforeEach } from 'vitest';
import { mockApiCall } from './helpers/mockApiCall';

const { serverApi } = await import('../src/api/api');

let getJson: ReturnType<typeof vi.fn>;
let postJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson, postJson } = mockApiCall());
});

describe('serverApi.init', () => {
  it('asks for no project when none is wanted', async () => {
    getJson.mockResolvedValueOnce({ ready: true, projects: ['A'] });

    expect(await serverApi.init()).toEqual({ ready: true, projects: ['A'] });
    expect(getJson).toHaveBeenCalledWith('/init', {}, { label: 'init' });
  });

  it('sends the wanted project', async () => {
    getJson.mockResolvedValueOnce({ ready: false });

    await serverApi.init('Alpha');

    expect(getJson).toHaveBeenCalledWith('/init', { project: 'Alpha' }, { label: 'init' });
  });

  it('returns null when the call failed', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await serverApi.init()).toBeNull();
  });
});

describe('serverApi.cors', () => {
  it('asks quietly and returns the origins', async () => {
    getJson.mockResolvedValueOnce({ origins: ['http://a'] });

    expect(await serverApi.cors()).toEqual({ origins: ['http://a'] });
    expect(getJson).toHaveBeenCalledWith('/cors', {}, { label: 'cors', silent: true });
  });

  it('returns null when the call failed', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await serverApi.cors()).toBeNull();
  });
});

describe('serverApi.ensureJupyter', () => {
  it('sends the folder and the theme, and returns the port', async () => {
    postJson.mockResolvedValueOnce({ port: 8888, root_dir: '/data' });

    expect(await serverApi.ensureJupyter('/data', true)).toEqual({ port: 8888, root_dir: '/data' });
    expect(postJson).toHaveBeenCalledWith('/jupyter/ensure', { root_dir: '/data', dark: true },
      { label: 'jupyter ensure' });
  });

  it('returns null when the call failed', async () => {
    postJson.mockResolvedValueOnce(null);

    expect(await serverApi.ensureJupyter('/data', false)).toBeNull();
  });
});
