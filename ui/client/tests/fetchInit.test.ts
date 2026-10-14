import { describe, it, expect, vi, beforeEach } from 'vitest';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetch = vi.fn();
global.fetch = mockFetch;

const { fetchInit } = await import('../src/io/fetchInit');
const { useProjectStore, NO_PROJECT } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');
const { useToolkitStore } = await import('../src/stores/useToolkitStore');
const { useServerConstants } = await import('../src/stores/useServerConstants');
const { useServerUserStore } = await import('../src/stores/useServerUserStore');

const jsonAnswer = (data: any) => ({
  ok: true,
  headers: { get: () => 'application/json' },
  json: () => Promise.resolve(data),
});

const fullAnswer = {
  ready: true,
  projects: [{ name: 'Alpha' }, { name: 'Beta' }],
  datatypes: { STRING: 'string' },
  toolkitDocs: [{ toolkit: 'LSM', desc: { classpath: 'lsm.cls', description: 'Linear' } }],
  project: { name: 'Alpha', documents: [] },
  server: { heraVersion: '9.9.9', username: 'eran', inDocker: true },
};

beforeEach(() => {
  vi.clearAllMocks();
  useProjectListStore.setState({ projectNames: [] });
  useToolkitStore.setState({ toolkits: [] });
  useServerConstants.setState({ dataTypes: {} });
  useServerUserStore.setState({ username: null, inDocker: false, heraVersion: null });
  useProjectStore.setState({ currProjectName: NO_PROJECT, currProject: null, editedDocs: {} });
});

describe('fetchInit', () => {
  it('fills every store from one call', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer(fullAnswer));

    const data = await fetchInit('Alpha');

    expect(mockFetch).toHaveBeenCalledTimes(1);
    expect(data?.ready).toBe(true);
    expect(useProjectListStore.getState().projectNames).toEqual([{ name: 'Alpha' }, { name: 'Beta' }]);
    expect(useServerConstants.getState().dataTypes).toEqual({ STRING: 'string' });
    expect(useToolkitStore.getState().toolkits).toEqual([
      { toolkit: 'LSM', cls: 'lsm.cls', description: 'Linear' },
    ]);
    expect(useServerUserStore.getState()).toMatchObject({
      username: 'eran', inDocker: true, heraVersion: '9.9.9',
    });
    expect(useProjectStore.getState().currProjectName).toBe('Alpha');
    expect(useProjectStore.getState().currProject).toEqual({ name: 'Alpha', documents: [] });
  });

  it('sends the wanted project name', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer(fullAnswer));

    await fetchInit('My Project');

    expect(mockFetch.mock.calls[0][0]).toBe('http://test/init?project=My%20Project');
  });

  it('asks for no project when the URL has none', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer(fullAnswer));

    await fetchInit(undefined);

    expect(mockFetch.mock.calls[0][0]).toBe('http://test/init');
  });

  it('fills nothing while the server is still warming', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({ ready: false }));

    const data = await fetchInit('Alpha');

    expect(data?.ready).toBe(false);
    expect(useProjectListStore.getState().projectNames).toEqual([]);
    expect(useProjectStore.getState().currProjectName).toBe(NO_PROJECT);
  });

  it('takes the pieces that came and skips the ones that failed', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer({
      ready: true,
      projects: [{ name: 'Alpha' }],
      project: { name: 'Alpha', documents: [] },
      problems: { toolkitDocs: 'KeyError: classpath' },
    }));

    const data = await fetchInit('Alpha');

    expect(data?.problems).toEqual({ toolkitDocs: 'KeyError: classpath' });
    expect(useToolkitStore.getState().toolkits).toEqual([]);
    expect(useProjectStore.getState().currProject).toEqual({ name: 'Alpha', documents: [] });
  });

  it('answers null when the server has no init route', async () => {
    mockFetch.mockResolvedValueOnce({
      ok: true,
      headers: { get: () => 'text/html' },
      json: () => Promise.resolve({}),
    });

    expect(await fetchInit('Alpha')).toBeNull();
    expect(useProjectListStore.getState().projectNames).toEqual([]);
  });
});
