/// <reference types="node" />
import { describe, it, expect, afterAll } from 'vitest';
import { SHARED_SERVER_URL } from './mockFactories';

const execPython = async (code: string) => {
  const r = await fetch(`${SHARED_SERVER_URL}/exec`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ code }),
  });
  const text = await r.text();
  if (!r.ok) {
    throw new Error(`exec failed (${r.status}): ${text}`);
  }
  return JSON.parse(text).data;
};

describe('Add Project API integration', () => {
  afterAll(async () => {
    await execPython(`
from hera.datalayer.project import deleteProject
try:
    deleteProject('APITestProject')
except Exception:
    pass
result = 'ok'
`);
  }, 15000);

  it('server is healthy', async () => {
    const r = await fetch(`${SHARED_SERVER_URL}/healthz`);
    expect(r.ok).toBe(true);
    const data = await r.json();
    expect(data.status).toBe('ok');
  });

  it('APITestProject does not exist initially', async () => {
    const data = await execPython(`
from hera.datalayer.project import getProjectList
result = [{"name": p} for p in getProjectList()]
    `);
    expect(data).not.toContainEqual({ name: 'APITestProject' });
  });

  it('created project appears in project list', async () => {
    const createResult = await execPython(`
import os
from types import SimpleNamespace
from hera.utils.data.CLI import project_create
from hera.datalayer.project import getProjectList, Project

project_create(SimpleNamespace(
  projectName='APITestProject',
  directory=os.path.join(os.getcwd(), 'projects', 'APITestProject'),
  loadRepositories=False,
  overwrite=False))

Project(projectName='APITestProject',
        filesDirectory=os.path.join(os.getcwd(), 'projects', 'APITestProject'))

result = [{"name": p} for p in getProjectList()]
    `);

    expect(createResult).toContainEqual({ name: 'APITestProject' });
  });

  it('project persists after simulated refresh (fresh getProjectList call)', async () => {
    const data = await execPython(`
from hera.datalayer.project import getProjectList
result = [{"name": p} for p in getProjectList()]
    `);

    expect(data).toContainEqual({ name: 'APITestProject' });
  });

  it('project has documents after refresh', async () => {
    const data = await execPython(`
from hera.datalayer import All
docs = All.getDocumentsAsDict('APITestProject', with_id=True)
result = docs['documents']
    `);

    expect(Array.isArray(data)).toBe(true);
    expect(data.length).toBeGreaterThan(0);
  });
}, 60000);
