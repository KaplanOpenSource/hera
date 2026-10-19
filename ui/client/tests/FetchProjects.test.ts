import { describe, it, expect, vi, beforeEach } from 'vitest';
import { resolveProjectFromUrl, fetchProjectsNames, fetchProjectDetails, fetchProjectData } from '../src/io/FetchProjects';
import { useProjectStore } from '../src/stores/useProjectStore';
import { useProjectListStore } from '../src/stores/useProjectListStore';
import { useToolkitStore } from '../src/stores/useToolkitStore';
import { mockApiCall } from './helpers/mockApiCall';

let getJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson } = mockApiCall());
});

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
    useProjectListStore.getState().setProjectNames([]);
    vi.mocked(getJson).mockReset();
  });

  it('asks the project.names entrypoint', async () => {
    vi.mocked(getJson).mockResolvedValue({ projects: [{ name: 'Alpha' }] });

    await fetchProjectsNames();

    expect(vi.mocked(getJson).mock.calls[0][0]).toBe('/project/names');
  });

  it('sets project names in store on success', async () => {
    vi.mocked(getJson).mockResolvedValue({ projects: [{ name: 'Alpha' }, { name: 'Beta' }] });

    await fetchProjectsNames();

    expect(useProjectListStore.getState().projectNames).toEqual([{ name: 'Alpha' }, { name: 'Beta' }]);
  });

  it('does not update store when the call failed', async () => {
    vi.mocked(getJson).mockResolvedValue(null);

    await fetchProjectsNames();

    expect(useProjectListStore.getState().projectNames).toEqual([]);
  });

  it('store not updated while response is in flight', async () => {
    let resolve!: (v: any) => void;
    vi.mocked(getJson).mockImplementationOnce(() =>
      new Promise(r => { resolve = r; })
    );

    const promise = fetchProjectsNames();
    expect(useProjectListStore.getState().projectNames).toEqual([]);

    resolve({ projects: [{ name: 'Gamma' }] });
    await promise;

    expect(useProjectListStore.getState().projectNames).toEqual([{ name: 'Gamma' }]);
  });
});

describe('fetchProjectDetails', () => {
  beforeEach(() => {
    useProjectStore.setState({
      currProjectName: 'TestProject',
      currProject: null,
    });
    vi.mocked(getJson).mockReset();
  });

  it('sets current project in store on success', async () => {
    const project = { name: 'TestProject', documents: [] };
    vi.mocked(getJson).mockResolvedValueOnce({ project });

    await fetchProjectDetails('TestProject');

    expect(useProjectStore.getState().currProject).toEqual(project);
  });

  it('does not update store when the call failed', async () => {
    vi.mocked(getJson).mockResolvedValueOnce(null);

    await fetchProjectDetails('TestProject');

    expect(useProjectStore.getState().currProject).toBeNull();
  });

  it('skips update if project name changed during fetch', async () => {
    const project = { name: 'TestProject', documents: [] };
    vi.mocked(getJson).mockImplementationOnce(async () => {
      useProjectStore.setState({ currProjectName: 'OtherProject' });
      return { project };
    });

    await fetchProjectDetails('TestProject');

    expect(useProjectStore.getState().currProject).toBeNull();
  });

  it('stale response discarded when user switches projects during fetch', async () => {
    let resolveAlpha!: (v: any) => void;
    let resolveBeta!: (v: any) => void;
    vi.mocked(getJson)
      .mockImplementationOnce(() => new Promise(r => { resolveAlpha = r; }))
      .mockImplementationOnce(() => new Promise(r => { resolveBeta = r; }));

    // User loads Alpha, then quickly switches to Beta
    const alphaPromise = fetchProjectDetails('Alpha');
    useProjectStore.setState({ currProjectName: 'Beta' });
    const betaPromise = fetchProjectDetails('Beta');

    // Beta responds first
    resolveBeta({ project: { name: 'Beta', documents: [{ x: 1 }] } });
    await betaPromise;
    expect(useProjectStore.getState().currProject).toEqual({ name: 'Beta', documents: [{ x: 1 }] });

    // Alpha responds late — should be discarded since currProjectName is now Beta
    resolveAlpha({ project: { name: 'Alpha', documents: [] } });
    await alphaPromise;
    expect(useProjectStore.getState().currProject).toEqual({ name: 'Beta', documents: [{ x: 1 }] });
  });

  it('sends the project name to the project.documents entrypoint', async () => {
    vi.mocked(getJson).mockResolvedValueOnce({ project: { name: 'TestProject', documents: [] } });

    await fetchProjectDetails('TestProject');

    const [path, params] = vi.mocked(getJson).mock.calls[0];
    expect(path).toBe('/project/documents');
    expect(params).toEqual({ projectName: 'TestProject' });
  });

  it('asks for no toasts when silent', async () => {
    vi.mocked(getJson).mockResolvedValueOnce({ project: { name: 'TestProject', documents: [] } });

    await fetchProjectDetails('TestProject', true);

    expect(vi.mocked(getJson).mock.calls[0][2].silent).toBe(true);
  });

  it('store not updated while response is in flight', async () => {
    const project = { name: 'TestProject', documents: [] };
    let resolve!: (v: any) => void;
    vi.mocked(getJson).mockImplementationOnce(() =>
      new Promise(r => { resolve = r; })
    );

    const promise = fetchProjectDetails('TestProject');
    expect(useProjectStore.getState().currProject).toBeNull();

    resolve({ project });
    await promise;

    expect(useProjectStore.getState().currProject).toEqual(project);
  });
});

describe('fetchProjectData', () => {
  beforeEach(() => {
    useProjectStore.setState({
      currProjectName: 'TestProject',
      currProject: null,
    });
    useToolkitStore.setState({ toolkits: [] });
    vi.mocked(getJson).mockReset();
  });

  const answerBoth = () => {
    vi.mocked(getJson).mockImplementation(async (path: string) => {
      if (path === '/project/toolkits') {
        return {
          toolkitDocs: [
            { toolkit: 'LSM', desc: { classpath: 'lsm.cls', description: 'Linear' } },
            { toolkit: 'GIS', desc: { classpath: 'gis.cls' } },
          ],
        };
      }
      return { project: { name: 'TestProject', documents: [] } };
    });
  };

  it('sets toolkits and project in store on success', async () => {
    answerBoth();

    await fetchProjectData('TestProject');

    expect(useToolkitStore.getState().toolkits).toEqual([
      { toolkit: 'LSM', cls: 'lsm.cls', description: 'Linear' },
      { toolkit: 'GIS', cls: 'gis.cls' },
    ]);
    expect(useProjectStore.getState().currProject).toEqual({ name: 'TestProject', documents: [] });
  });

  it('asks both entrypoints', async () => {
    answerBoth();

    await fetchProjectData('TestProject');

    const paths = vi.mocked(getJson).mock.calls.map(c => c[0]);
    expect(paths).toEqual(['/project/toolkits', '/project/documents']);
  });

  it('does not update store when both calls failed', async () => {
    vi.mocked(getJson).mockResolvedValue(null);

    await fetchProjectData('TestProject');

    expect(useToolkitStore.getState().toolkits).toEqual([]);
    expect(useProjectStore.getState().currProject).toBeNull();
  });

  it('keeps the toolkits when only the project call failed', async () => {
    vi.mocked(getJson).mockImplementation(async (path: string) => {
      if (path === '/project/toolkits') {
        return { toolkitDocs: [{ toolkit: 'T1', desc: { classpath: 't1.cls' } }] };
      }
      return null;
    });

    await fetchProjectData('TestProject');

    expect(useToolkitStore.getState().toolkits).toEqual([{ toolkit: 'T1', cls: 't1.cls' }]);
    expect(useProjectStore.getState().currProject).toBeNull();
  });

  it('store not updated while response is in flight', async () => {
    let resolve!: (v: any) => void;
    vi.mocked(getJson)
      .mockImplementationOnce(async () => ({ toolkitDocs: [{ toolkit: 'T1', desc: { classpath: 't1.cls' } }] }))
      .mockImplementationOnce(() => new Promise(r => { resolve = r; }));

    const promise = fetchProjectData('TestProject');
    expect(useToolkitStore.getState().toolkits).toEqual([]);
    expect(useProjectStore.getState().currProject).toBeNull();

    resolve({ project: { name: 'TestProject', documents: [] } });
    await promise;

    expect(useToolkitStore.getState().toolkits).toEqual([{ toolkit: 'T1', cls: 't1.cls' }]);
    expect(useProjectStore.getState().currProject).toEqual({ name: 'TestProject', documents: [] });
  });
});
