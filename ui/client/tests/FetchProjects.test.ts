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

  it('works with delayed response', async () => {
    vi.mocked(fetchPython).mockImplementationOnce(() =>
      new Promise(resolve => setTimeout(() => resolve({
        data: { projects: [{ name: 'Gamma' }] },
        problem: undefined,
      }), 50))
    );

    await fetchProjectsNames();

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

  it('works with delayed response', async () => {
    const project = { name: 'TestProject', documents: [] };
    vi.mocked(execPython).mockImplementationOnce(() =>
      new Promise(resolve => setTimeout(() => resolve({
        data: JSON.stringify(project),
        problem: undefined,
      }), 50))
    );

    await fetchProjectDetails('TestProject');

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

  it('works with delayed response', async () => {
    vi.mocked(execPython).mockImplementationOnce(() =>
      new Promise(resolve => setTimeout(() => resolve({
        data: [{ toolkit: 'T1', desc: { classpath: 't1.cls' } }],
        problem: undefined,
      }), 50))
    );

    await fetchToolkits('TestProject');

    expect(useProjectStore.getState().toolkits).toEqual([
      { toolkit: 'T1', cls: 't1.cls' },
    ]);
  });
});
