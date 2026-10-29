import { describe, it, expect, vi, beforeEach } from 'vitest';
import { mockApiCall } from './helpers/mockApiCall';
import { ApiCall } from '../src/api/ApiCall';

const { fetchInit } = await import('../src/io/fetchInit');
const { useProjectStore, NO_PROJECT } = await import('../src/stores/useProjectStore');
const { useProjectListStore } = await import('../src/stores/useProjectListStore');
const { useToolkitStore } = await import('../src/stores/useToolkitStore');
const { useServerConstants } = await import('../src/stores/useServerConstants');
const { useServerUserStore } = await import('../src/stores/useServerUserStore');

const fullAnswer = {
  ready: true,
  projects: [{ name: 'Alpha' }, { name: 'Beta' }],
  datatypes: { STRING: 'string' },
  toolkitDocs: [{ toolkit: 'LSM', desc: { classpath: 'lsm.cls', description: 'Linear' } }],
  project: { name: 'Alpha', documents: [] },
  server: { heraVersion: '9.9.9', username: 'eran', inDocker: true },
};

let getJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson } = mockApiCall());
  useProjectListStore.setState({ projectNames: [] });
  useToolkitStore.setState({ toolkits: [] });
  useServerConstants.setState({ dataTypes: {} });
  useServerUserStore.setState({ username: null, inDocker: false, heraVersion: null });
  useProjectStore.setState({ currProjectName: NO_PROJECT, currProject: null, editedDocs: {} });
});

describe('fetchInit', () => {
  it('fills every store from one call', async () => {
    getJson.mockResolvedValueOnce(fullAnswer);

    const data = await fetchInit('Alpha');

    expect(getJson).toHaveBeenCalledTimes(1);
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

  it('starts notebook detection for the project it loaded', async () => {
    const postJson = vi.spyOn(ApiCall.prototype as any, 'postJson').mockResolvedValue(null);
    const configDoc = { _id: { $oid: 'cfg1' }, type: 'Alpha__config__', desc: { filesDirectory: '/data/a' } };
    getJson.mockResolvedValueOnce({ ...fullAnswer, project: { name: 'Alpha', documents: [configDoc] } });

    await fetchInit('Alpha');

    expect(postJson).toHaveBeenCalledWith('/project/detectNotebooks',
      { projectName: 'Alpha', filesDirectory: '/data/a' },
      { label: 'detect notebooks' });
  });

  it('sends the wanted project name', async () => {
    getJson.mockResolvedValueOnce(fullAnswer);

    await fetchInit('My Project');

    expect(getJson).toHaveBeenCalledWith('/init', { project: 'My Project' }, { label: 'init' });
  });

  it('asks for no project when the URL has none', async () => {
    getJson.mockResolvedValueOnce(fullAnswer);

    await fetchInit(undefined);

    expect(getJson).toHaveBeenCalledWith('/init', {}, { label: 'init' });
  });

  it('fills nothing while the server is still warming', async () => {
    getJson.mockResolvedValueOnce({ ready: false });

    const data = await fetchInit('Alpha');

    expect(data?.ready).toBe(false);
    expect(useProjectListStore.getState().projectNames).toEqual([]);
    expect(useProjectStore.getState().currProjectName).toBe(NO_PROJECT);
  });

  it('takes the pieces that came and skips the ones that failed', async () => {
    getJson.mockResolvedValueOnce({
      ready: true,
      projects: [{ name: 'Alpha' }],
      project: { name: 'Alpha', documents: [] },
      problems: { toolkitDocs: 'KeyError: classpath' },
    });

    const data = await fetchInit('Alpha');

    expect(data?.problems).toEqual({ toolkitDocs: 'KeyError: classpath' });
    expect(useToolkitStore.getState().toolkits).toEqual([]);
    expect(useProjectStore.getState().currProject).toEqual({ name: 'Alpha', documents: [] });
  });

  it('answers null when the server has no init route', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await fetchInit('Alpha')).toBeNull();
    expect(useProjectListStore.getState().projectNames).toEqual([]);
  });
});
