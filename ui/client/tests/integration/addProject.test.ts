/// <reference types="node" />
import { describe, it, expect, beforeAll, afterAll } from 'vitest';
import { execSync } from 'child_process';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';

const TEST_SERVER_PORT = 8001;
const SERVER_URL = `http://localhost:${TEST_SERVER_PORT}`;
const TEST_DB_NAME = 'hera_test';
const DOCKER_NETWORK = 'hera-test-net';
const MONGO_CONTAINER = 'hera-test-mongo';
const SERVER_CONTAINER = 'hera-test-server';
const PROJECT_ROOT = path.resolve(__dirname, '../../../../');

let tmpDir: string;

const docker = (cmd: string) => execSync(`docker ${cmd}`, { encoding: 'utf-8' }).trim();

const execPython = async (code: string) => {
  const r = await fetch(`${SERVER_URL}/exec`, {
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

const waitForServer = async (timeoutMs = 20000) => {
  const start = Date.now();
  while (Date.now() - start < timeoutMs) {
    try {
      const r = await fetch(`${SERVER_URL}/healthz`);
      if (r.ok) return;
    } catch { /* not ready yet */ }
    await new Promise(r => setTimeout(r, 300));
  }
  throw new Error(`Server not ready after ${timeoutMs}ms`);
};

const cleanup = () => {
  for (const name of [SERVER_CONTAINER, MONGO_CONTAINER]) {
    try { docker(`rm -f ${name}`); } catch { /* ignore */ }
  }
  try { docker(`network rm ${DOCKER_NETWORK}`); } catch { /* ignore */ }
  if (tmpDir) {
    try { fs.rmSync(tmpDir, { recursive: true, force: true }); } catch { /* ignore */ }
  }
};

describe('Add Project integration', () => {
  beforeAll(async () => {
    // Clean up any leftovers from a previous failed run
    cleanup();

    tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'hera-integ-'));
    const pyheraDir = path.join(tmpDir, '.pyhera');
    fs.mkdirSync(pyheraDir);

    // Write pyhera config pointing to the test MongoDB container
    fs.writeFileSync(path.join(pyheraDir, 'config.json'), JSON.stringify({
      root: {
        dbName: TEST_DB_NAME,
        dbIP: MONGO_CONTAINER,
        username: 'hera',
        password: 'heracles',
      },
    }, null, 2));

    // Create Docker network
    docker(`network create ${DOCKER_NETWORK}`);

    // Start test MongoDB
    const mongoInitDir = path.join(PROJECT_ROOT, 'mongo-init.d');
    docker([
      `run -d --name ${MONGO_CONTAINER}`,
      `--network ${DOCKER_NETWORK}`,
      `-v ${mongoInitDir}:/docker-entrypoint-initdb.d`,
      `mongo:5.0`,
    ].join(' '));

    // Wait for MongoDB to be ready
    const mongoStart = Date.now();
    while (Date.now() - mongoStart < 15000) {
      try {
        const out = docker(`exec ${MONGO_CONTAINER} mongosh --eval "db.runCommand({ping:1})" --quiet`);
        if (out.includes('ok')) break;
      } catch { /* not ready */ }
      await new Promise(r => setTimeout(r, 300));
    }

    // Start test server from the hera-server image
    docker([
      `run -d --name ${SERVER_CONTAINER}`,
      `--network ${DOCKER_NETWORK}`,
      `-p 127.0.0.1:${TEST_SERVER_PORT}:8000`,
      `-v ${PROJECT_ROOT}:/app`,
      `-v ${pyheraDir}:/root/.pyhera`,
      `hera-server`,
      `python ui/server/server.py --cors`,
    ].join(' '));

    await waitForServer();
  }, 30000);

  afterAll(() => {
    cleanup();
  }, 15000);

  it('server is healthy', async () => {
    const r = await fetch(`${SERVER_URL}/healthz`);
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
