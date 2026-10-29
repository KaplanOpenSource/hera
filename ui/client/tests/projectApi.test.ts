import { describe, it, expect, vi, beforeEach } from 'vitest';
import { mockApiCall } from './helpers/mockApiCall';

const { projectApi } = await import('../src/api/api');

let getJson: ReturnType<typeof vi.fn>;
let postJson: ReturnType<typeof vi.fn>;

beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson, postJson } = mockApiCall());
});

describe('projectApi.names', () => {
  it('asks for the names and returns the answer', async () => {
    getJson.mockResolvedValueOnce({ projects: ['A', 'B'] });

    expect(await projectApi.names()).toEqual({ projects: ['A', 'B'] });
    expect(getJson).toHaveBeenCalledWith('/project/names', {}, { label: 'project names' });
  });

  it('returns null when the call failed', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await projectApi.names()).toBeNull();
  });
});

describe('projectApi.documents', () => {
  it('sends the project and the label, not silent by default', async () => {
    getJson.mockResolvedValueOnce({ project: { projectName: 'A' } });

    expect(await projectApi.documents('A', 'load A')).toEqual({ project: { projectName: 'A' } });
    expect(getJson).toHaveBeenCalledWith('/project/documents', { projectName: 'A' },
      { label: 'load A', silent: false });
  });

  it('can ask quietly', async () => {
    getJson.mockResolvedValueOnce({});

    await projectApi.documents('A', 'reload', true);

    expect(getJson).toHaveBeenCalledWith('/project/documents', { projectName: 'A' },
      { label: 'reload', silent: true });
  });

  it('returns null when the call failed', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await projectApi.documents('A', 'load A')).toBeNull();
  });
});

describe('projectApi.toolkits', () => {
  it('asks for the toolkits and returns the answer', async () => {
    getJson.mockResolvedValueOnce({ toolkitDocs: [] });

    expect(await projectApi.toolkits()).toEqual({ toolkitDocs: [] });
    expect(getJson).toHaveBeenCalledWith('/project/toolkits', {}, { label: 'toolkits' });
  });

  it('returns null when the call failed', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await projectApi.toolkits()).toBeNull();
  });
});

describe('projectApi.create', () => {
  it('sends the name, the folder and the repositories flag', async () => {
    postJson.mockResolvedValueOnce({ project: { projectName: 'A' } });

    expect(await projectApi.create('A', '/data/a', true)).toEqual({ project: { projectName: 'A' } });
    expect(postJson).toHaveBeenCalledWith('/project/create',
      { projectName: 'A', filesDirectory: '/data/a', loadRepositories: true },
      { label: 'create project A' });
  });

  it('returns null when the call failed', async () => {
    postJson.mockResolvedValueOnce(null);

    expect(await projectApi.create('A', '', false)).toBeNull();
  });
});

describe('projectApi.detectNotebooks', () => {
  it('sends the name and the folder', async () => {
    postJson.mockResolvedValueOnce({ added: 2 });

    expect(await projectApi.detectNotebooks('A', '/data/a')).toEqual({ added: 2 });
    expect(postJson).toHaveBeenCalledWith('/project/detectNotebooks',
      { projectName: 'A', filesDirectory: '/data/a' },
      { label: 'detect notebooks' });
  });

  it('returns null when the call failed', async () => {
    postJson.mockResolvedValueOnce(null);

    expect(await projectApi.detectNotebooks('A', '/data/a')).toBeNull();
  });
});

describe('projectApi.delete', () => {
  it('sends the name and the delete files flag', async () => {
    postJson.mockResolvedValueOnce({ projects: [] });

    expect(await projectApi.delete('A', true)).toEqual({ projects: [] });
    expect(postJson).toHaveBeenCalledWith('/project/delete',
      { projectName: 'A', deleteFiles: true },
      { label: 'delete project A' });
  });

  it('returns null when the call failed', async () => {
    postJson.mockResolvedValueOnce(null);

    expect(await projectApi.delete('A', false)).toBeNull();
  });
});
