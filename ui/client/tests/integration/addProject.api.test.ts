/// <reference types="node" />
import { describe, it, expect, beforeAll, afterAll } from 'vitest';
import { SHARED_SERVER_URL } from './mockFactories';

const PROJECT = 'APITestProject';

const getJson = async (path: string) => {
  const r = await fetch(`${SHARED_SERVER_URL}${path}`);
  const text = await r.text();
  if (!r.ok) {
    throw new Error(`GET ${path} failed (${r.status}): ${text}`);
  }
  return JSON.parse(text);
};

const postJson = async (path: string, body: object) => {
  const r = await fetch(`${SHARED_SERVER_URL}${path}`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });
  const text = await r.text();
  if (!r.ok) {
    throw new Error(`POST ${path} failed (${r.status}): ${text}`);
  }
  return JSON.parse(text);
};

const projectNames = async (): Promise<{ name: string }[]> => {
  return (await getJson('/project/names')).projects;
};

// The project routes answer 503 until hera has warmed up.
const waitUntilReady = async () => {
  const deadline = Date.now() + 30_000;
  while (Date.now() < deadline) {
    if ((await getJson('/init')).ready) {
      return;
    }
    await new Promise(r => setTimeout(r, 500));
  }
  throw new Error('the server did not warm up within 30s');
};

describe('Add Project API integration', () => {
  beforeAll(waitUntilReady, 35000);

  afterAll(async () => {
    if ((await projectNames()).some(p => p.name === PROJECT)) {
      await postJson('/project/delete', { projectName: PROJECT, deleteFiles: true });
    }
  }, 15000);

  it('server is healthy', async () => {
    const data = await getJson('/healthz');
    expect(data.status).toBe('ok');
  });

  it('APITestProject does not exist initially', async () => {
    expect(await projectNames()).not.toContainEqual({ name: PROJECT });
  });

  it('created project appears in project list', async () => {
    const created = await postJson('/project/create', { projectName: PROJECT, loadRepositories: false });

    expect(created.projects).toContainEqual({ name: PROJECT });
  });

  it('project persists after simulated refresh (fresh project list call)', async () => {
    expect(await projectNames()).toContainEqual({ name: PROJECT });
  });

  it('project has documents after refresh', async () => {
    const { project } = await getJson(`/project/documents?projectName=${PROJECT}`);

    expect(Array.isArray(project.documents)).toBe(true);
    expect(project.documents.length).toBeGreaterThan(0);
  });
}, 60000);
