import { describe, it, expect, vi, beforeEach } from 'vitest';
import { resolveProjectFromUrl, fetchProjectsNames } from '../src/io/FetchProjects';
import { useProjectStore } from '../src/stores/useProjectStore';

vi.mock('../src/io/fetchPython', () => ({
  fetchPython: vi.fn(),
}));

import { fetchPython } from '../src/io/fetchPython';

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
});
