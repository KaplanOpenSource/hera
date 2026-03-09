import { describe, it, expect, vi, beforeEach } from 'vitest';
import { resolveProjectFromUrl, fetchProjectsNames, fetchProjectDetails, fetchToolkits } from '../src/io/FetchProjects';
import { useProjectStore, NO_PROJECT } from '../src/stores/useProjectStore';

vi.mock('../src/io/fetchPython', () => ({
  fetchPython: vi.fn(),
}));

import { fetchPython } from '../src/io/fetchPython';

vi.mock('../src/io/execPython', () => ({
  execPython: vi.fn(),
}));

import { execPython } from '../src/io/execPython';

const projects = [
  { name: 'Alpha' },
  { name: 'Beta' },
  { name: 'My Project' },
];

describe('resolveProjectFromUrl', () => {
  it('returns matching project name for a valid encoded URL param', () => {
    expect(resolveProjectFromUrl(encodeURIComponent('My Project'), projects)).toBe('My Project');
  });

  it('returns undefined when URL param does not match any project', () => {
    expect(resolveProjectFromUrl('NonExistent', projects)).toBeUndefined();
  });
});

describe('fetchProjectsNames', () => {
  beforeEach(() => {
    useProjectStore.getState().setProjectNames([]);
    vi.mocked(fetchPython).mockReset();
  });

  it('sets project names in store on success', async () => {
    vi.mocked(fetchPython).mockResolvedValue({
      data: { projects: [{ name: 'Alpha' }, { name: 'Beta' }] },
      problem: undefined,
    });

    await fetchProjectsNames();

    const names = useProjectStore.getState().projectNames;
    expect(names).toEqual([{ name: 'Alpha' }, { name: 'Beta' }]);
  });

  it('does not update store when there is a problem', async () => {
    vi.mocked(fetchPython).mockResolvedValue({
      data: undefined,
      problem: 'connection failed',
    });

    await fetchProjectsNames();

    const names = useProjectStore.getState().projectNames;
    expect(names).toEqual([]);
  });

  it('store not updated while response is in flight', async () => {
    let resolve!: (v: any) => void;
    vi.mocked(fetchPython).mockImplementationOnce(() =>
      new Promise(r => { resolve = r; })
    );

    const promise = fetchProjectsNames();
    expect(useProjectStore.getState().projectNames).toEqual([]);

    resolve({ data: { projects: [{ name: 'Gamma' }] }, problem: undefined });
    await promise;

    expect(useProjectStore.getState().projectNames).toEqual([{ name: 'Gamma' }]);
  });
});

describe('fetchProjectDetails', () => {
  beforeEach(() => {
    useProjectStore.setState({
      currProjectName: 'TestProject',
      currProject: null,
    });
    vi.mocked(execPython).mockReset();
  });

  it('sets current project in store on success', async () => {
    const project = { name: 'TestProject', documents: [] };
    vi.mocked(execPython).mockResolvedValueOnce({
      data: JSON.stringify(project),
      problem: undefined,
    });

    await fetchProjectDetails('TestProject');

    expect(useProjectStore.getState().currProject).toEqual(project);
  });

  it('does not update store when there is a problem', async () => {
    vi.mocked(execPython).mockResolvedValueOnce({
      data: undefined,
      problem: 'error',
    });

    await fetchProjectDetails('TestProject');

    expect(useProjectStore.getState().currProject).toBeNull();
  });

  it('skips update if project name changed during fetch', async () => {
    const project = { name: 'TestProject', documents: [] };
    vi.mocked(execPython).mockImplementationOnce(async () => {
      useProjectStore.setState({ currProjectName: 'OtherProject' });
      return { data: JSON.stringify(project), problem: undefined };
    });

    await fetchProjectDetails('TestProject');

    expect(useProjectStore.getState().currProject).toBeNull();
  });

  it('stale response discarded when user switches projects during fetch', async () => {
    let resolveAlpha!: (v: any) => void;
    let resolveBeta!: (v: any) => void;
    vi.mocked(execPython)
      .mockImplementationOnce(() => new Promise(r => { resolveAlpha = r; }))
      .mockImplementationOnce(() => new Promise(r => { resolveBeta = r; }));

    // User loads Alpha, then quickly switches to Beta
    const alphaPromise = fetchProjectDetails('Alpha');
    useProjectStore.setState({ currProjectName: 'Beta' });
    const betaPromise = fetchProjectDetails('Beta');

    // Beta responds first
    resolveBeta({ data: JSON.stringify({ name: 'Beta', documents: [{ x: 1 }] }), problem: undefined });
    await betaPromise;
    expect(useProjectStore.getState().currProject).toEqual({ name: 'Beta', documents: [{ x: 1 }] });

    // Alpha responds late — should be discarded since currProjectName is now Beta
    resolveAlpha({ data: JSON.stringify({ name: 'Alpha', documents: [] }), problem: undefined });
    await alphaPromise;
    expect(useProjectStore.getState().currProject).toEqual({ name: 'Beta', documents: [{ x: 1 }] });
  });

  it('sends correct Python code', async () => {
    vi.mocked(execPython).mockResolvedValueOnce({
      data: JSON.stringify({ name: 'TestProject', documents: [] }),
      problem: undefined,
    });

    await fetchProjectDetails('TestProject');

    const code = vi.mocked(execPython).mock.calls[0][0];
    expect(code).toContain("All.getDocumentsAsDict('TestProject'");
    expect(code).toContain('json.dumps');
  });

  it('store not updated while response is in flight', async () => {
    const project = { name: 'TestProject', documents: [] };
    let resolve!: (v: any) => void;
    vi.mocked(execPython).mockImplementationOnce(() =>
      new Promise(r => { resolve = r; })
    );

    const promise = fetchProjectDetails('TestProject');
    expect(useProjectStore.getState().currProject).toBeNull();

    resolve({ data: JSON.stringify(project), problem: undefined });
    await promise;

    expect(useProjectStore.getState().currProject).toEqual(project);
  });
});

describe('fetchToolkits', () => {
  beforeEach(() => {
    useProjectStore.setState({ toolkits: [] });
    vi.mocked(execPython).mockReset();
  });

  it('sets toolkits in store on success', async () => {
    vi.mocked(execPython).mockResolvedValueOnce({
      data: [
        { toolkit: 'LSM', desc: { classpath: 'lsm.cls', description: 'Linear' } },
        { toolkit: 'GIS', desc: { classpath: 'gis.cls' } },
      ],
      problem: undefined,
    });

    await fetchToolkits('TestProject');

    const toolkits = useProjectStore.getState().toolkits;
    expect(toolkits).toEqual([
      { toolkit: 'LSM', cls: 'lsm.cls', description: 'Linear' },
      { toolkit: 'GIS', cls: 'gis.cls' },
    ]);
  });

  it('does not update store when there is a problem', async () => {
    vi.mocked(execPython).mockResolvedValueOnce({
      data: undefined,
      problem: 'error',
    });

    await fetchToolkits('TestProject');

    expect(useProjectStore.getState().toolkits).toEqual([]);
  });

  it('store not updated while response is in flight', async () => {
    let resolve!: (v: any) => void;
    vi.mocked(execPython).mockImplementationOnce(() =>
      new Promise(r => { resolve = r; })
    );

    const promise = fetchToolkits('TestProject');
    expect(useProjectStore.getState().toolkits).toEqual([]);

    resolve({ data: [{ toolkit: 'T1', desc: { classpath: 't1.cls' } }], problem: undefined });
    await promise;

    expect(useProjectStore.getState().toolkits).toEqual([
      { toolkit: 'T1', cls: 't1.cls' },
    ]);
  });
});
