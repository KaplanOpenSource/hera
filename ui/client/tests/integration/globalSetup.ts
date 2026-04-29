/// <reference types="node" />
import { execSync, spawn } from 'child_process';
import * as net from 'net';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import { checkDocker, startDockerEnv, stopDockerEnv, type DockerEnv } from './dockerSetup';

const PORT_FILE = path.join(os.tmpdir(), 'hera-integ-port');
const CLEANUP_FILE = path.join(os.tmpdir(), 'hera-integ-cleanup.json');
const findProjectRoot = (): string => {
  let dir = __dirname;
  while (dir !== path.dirname(dir)) {
    if (fs.existsSync(path.join(dir, 'ui/server/server.py'))) return dir;
    dir = path.dirname(dir);
  }
  if (fs.existsSync('/app/ui/server/server.py')) return '/app';
  throw new Error('Cannot find project root (looking for ui/server/server.py)');
};

const PROJECT_ROOT = findProjectRoot();

type CleanupManifest = {
  serverPid?: number;
  tmpDir?: string;
} & Partial<DockerEnv>;

// ---- Probing ----

const probeServer = async (port: number): Promise<boolean> => {
  try {
    const r = await fetch(`http://localhost:${port}/healthz`);
    if (!r.ok) return false;
    const data = await r.json();
    return data.status === 'ok';
  } catch {
    return false;
  }
};

const probeTcp = (host: string, port: number): Promise<boolean> =>
  new Promise(resolve => {
    const socket = new net.Socket();
    socket.setTimeout(2000);
    socket.on('connect', () => { socket.destroy(); resolve(true); });
    socket.on('error', () => { socket.destroy(); resolve(false); });
    socket.on('timeout', () => { socket.destroy(); resolve(false); });
    socket.connect(port, host);
  });

// ---- Start / stop: server as a local process ----

const stopProcess = (pid: number) => {
  try { process.kill(pid, 'SIGTERM'); } catch { /* already exited */ }
};

const findPython = (): string | null => {
  for (const cmd of ['python3', 'python']) {
    try {
      return execSync(`which ${cmd}`, { encoding: 'utf-8', stdio: 'pipe' }).trim() || null;
    } catch { /* not in PATH */ }
  }
  return null;
};

const findMongo = async (): Promise<number | null> => {
  for (const port of [27018, 27017]) {
    if (await probeTcp('127.0.0.1', port)) return port;
  }
  return null;
};

const startServerProcess = async (mongoPort: number) => {
  const python = findPython();
  if (!python) throw new Error('No Python binary found');

  const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'hera-integ-'));
  const pyheraDir = path.join(tmpDir, '.pyhera');
  fs.mkdirSync(pyheraDir);

  const dbIP = mongoPort === 27017 ? '127.0.0.1' : `127.0.0.1:${mongoPort}`;

  fs.writeFileSync(path.join(pyheraDir, 'config.json'), JSON.stringify({
    [os.userInfo().username]: {
      dbIP,
      dbName: 'olymp',
      username: 'hera',
      password: 'heracles',
    },
  }, null, 2));

  const proc = spawn(python, [
    'ui/server/server.py', '--cors', 'all', '-y', '--jupyter-port', '0',
  ], {
    cwd: PROJECT_ROOT,
    env: { ...process.env, HOME: tmpDir },
    stdio: ['pipe', 'pipe', 'pipe'],
  });

  const pid = proc.pid;
  if (!pid) {
    fs.rmSync(tmpDir, { recursive: true, force: true });
    throw new Error(`Failed to spawn server process with ${python}`);
  }

  const serverLog = path.join(tmpDir, 'server.log');
  const logStream = fs.createWriteStream(serverLog);
  proc.stdout?.pipe(logStream);
  proc.stderr?.pipe(logStream);

  let earlyExit: number | null = null;
  let stderr = '';
  proc.stderr?.on('data', (chunk: Buffer) => { stderr += chunk.toString(); });
  proc.on('exit', (code) => { earlyExit = code; });

  const port = 8000;
  const deadline = Date.now() + 30_000;
  while (Date.now() < deadline) {
    if (earlyExit !== null) {
      fs.rmSync(tmpDir, { recursive: true, force: true });
      throw new Error(
        `Server process exited early with code ${earlyExit}`
        + (stderr ? `\n${stderr.slice(-500)}` : ''),
      );
    }
    if (await probeServer(port)) break;
    await new Promise(r => setTimeout(r, 500));
  }
  if (!(await probeServer(port))) {
    stopProcess(pid);
    fs.rmSync(tmpDir, { recursive: true, force: true });
    throw new Error('Server process did not become ready within 30 s');
  }

  return { pid, port, tmpDir };
};

const stopServerProcess = (pid: number, tmpDir: string) => {
  stopProcess(pid);
  try { fs.rmSync(tmpDir, { recursive: true, force: true }); } catch { /* best effort */ }
};

// ---- Cleanup manifest (crash recovery) ----

const writeManifest = (m: CleanupManifest) =>
  fs.writeFileSync(CLEANUP_FILE, JSON.stringify(m));

const cleanupPreviousRun = () => {
  let m: CleanupManifest;
  try { m = JSON.parse(fs.readFileSync(CLEANUP_FILE, 'utf-8')); } catch { return; }

  if (m.serverPid) stopProcess(m.serverPid);
  if (m.serverContainerId || m.mongoContainerId) {
    stopDockerEnv(m as DockerEnv);
  }
  if (m.tmpDir) {
    try { fs.rmSync(m.tmpDir, { recursive: true, force: true }); } catch { /* best effort */ }
  }
  try { fs.unlinkSync(CLEANUP_FILE); } catch { /* already gone */ }
};

// ---- Signal-safe teardown registration ----

const registerCleanup = (stopFn: () => void) => {
  let cleaned = false;
  const cleanup = () => {
    if (cleaned) return;
    cleaned = true;
    stopFn();
    try { fs.unlinkSync(CLEANUP_FILE); } catch { /* ok */ }
    try { fs.unlinkSync(PORT_FILE); } catch { /* ok */ }
  };

  const onSigint = () => { cleanup(); process.exit(130); };
  const onSigterm = () => { cleanup(); process.exit(143); };
  process.on('SIGINT', onSigint);
  process.on('SIGTERM', onSigterm);
  process.on('exit', cleanup);

  return function teardown() {
    cleanup();
    process.removeListener('SIGINT', onSigint);
    process.removeListener('SIGTERM', onSigterm);
    process.removeListener('exit', cleanup);
  };
};

// ---- Orchestrator ----

export default async function setup() {
  cleanupPreviousRun();

  // Strategy 1: already-running server
  for (const port of [8001, 8000]) {
    if (await probeServer(port)) {
      console.log(`[globalSetup] Using existing server on port ${port}`);
      fs.writeFileSync(PORT_FILE, String(port));
      return function teardown() {
        try { fs.unlinkSync(PORT_FILE); } catch { /* ok */ }
      };
    }
  }

  // Strategy 2: mongo running + python available → start server as process
  const mongoPort = await findMongo();
  if (mongoPort && findPython()) {
    console.log(`[globalSetup] Starting server process (mongo on port ${mongoPort})`);
    const server = await startServerProcess(mongoPort);
    fs.writeFileSync(PORT_FILE, String(server.port));
    writeManifest({ serverPid: server.pid, tmpDir: server.tmpDir });
    return registerCleanup(() => stopServerProcess(server.pid, server.tmpDir));
  }

  // Strategy 3: Docker available → start mongo + server containers
  if (checkDocker()) {
    console.log('[globalSetup] Starting Docker environment on port 8001');
    const env = await startDockerEnv({
      network: `hera-test-net-${Date.now()}`,
      mongoContainer: `hera-test-mongo-${Date.now()}`,
      serverContainer: `hera-test-server-${Date.now()}`,
      serverPort: 8001,
      dbName: 'hera_test',
    });
    fs.writeFileSync(PORT_FILE, String(env.port));
    writeManifest(env);
    return registerCleanup(() => stopDockerEnv(env));
  }

  throw new Error(
    'No integration test environment available.\n'
    + 'Either:\n'
    + '  1. Start a hera server on port 8000 or 8001\n'
    + '  2. Have MongoDB running + Python with hera installed\n'
    + '  3. Install Docker with the hera-server image built',
  );
}
