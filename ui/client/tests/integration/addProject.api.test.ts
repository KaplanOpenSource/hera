/// <reference types="node" />
import { describe, it, expect, beforeAll, afterAll } from 'vitest';
import { startDockerEnv, type DockerEnv } from './dockerSetup';

let env: DockerEnv;

const execPython = async (code: string) => {
  const r = await fetch(`${env.serverUrl}/exec`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ code }),
  });
  const text = await r.text();
  if (!r.ok) {
    throw new Error(`exec failed (${r.status}): ${text}`);
  }
  return JSON.parse(text);
};

describe('Add Project API integration', () => {
  beforeAll(async () => {
    env = await startDockerEnv({
      network: 'hera-test-api-net',
      mongoContainer: 'hera-test-api-mongo',
      serverContainer: 'hera-test-api-server',
      serverPort: 8001,
      dbName: 'hera_test_api',
    });
  }, 30000);

  afterAll(() => {
    env?.cleanup();
  }, 15000);

  it('server is healthy', async () => {
    const r = await fetch(`${env.serverUrl}/healthz`);
    expect(r.ok).toBe(true);
    const data = await r.json();
    expect(data.status).toBe('ok');
  });

  it('project list is initially empty', async () => {
    const data = await execPython(`
from hera.datalayer.project import getProjectList
result = getProjectList()
    `);
    expect(data).toEqual([]);
  });

  it('created project appears in project list', async () => {
    const createResult = await execPython(`
import os
from types import SimpleNamespace
from hera.utils.data.CLI import project_create
from hera.datalayer.project import getProjectList, Project

project_create(SimpleNamespace(
  projectName='TestProject',
  directory=os.path.join(os.getcwd(), 'projects', 'TestProject'),
  loadRepositories=False,
  overwrite=False))

# Creates a config document in MongoDB so the project appears in getProjectList()
Project(projectName='TestProject',
        filesDirectory=os.path.join(os.getcwd(), 'projects', 'TestProject'))

result = [{"name": p} for p in getProjectList()]
    `);

    expect(createResult).toContainEqual({ name: 'TestProject' });
  });

  it('project persists after simulated refresh (fresh getProjectList call)', async () => {
    const data = await execPython(`
from hera.datalayer.project import getProjectList
result = [{"name": p} for p in getProjectList()]
    `);

    expect(data).toContainEqual({ name: 'TestProject' });
  });

  it('project has documents after refresh', async () => {
    const data = await execPython(`
from hera.datalayer import All
docs = All.getDocumentsAsDict('TestProject', with_id=True)
result = docs['documents']
    `);

    expect(Array.isArray(data)).toBe(true);
    expect(data.length).toBeGreaterThan(0);
  });
}, 60000);
