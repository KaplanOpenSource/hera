/// <reference types="node" />
import { spawn } from 'child_process';
import * as fs from 'fs';
import * as net from 'net';
import * as os from 'os';
import * as path from 'path';

const PORT_FILE = path.join(os.tmpdir(), 'hera-integ-port');
const MONGO_PORT = Number(process.env.HERA_INTEG_MONGO_PORT) || 27018;

// Ask the OS for a port nobody is using. Picking a fixed port is what let an
// earlier version talk to whatever server already happened to be listening.
const findFreePort = (): Promise<number> => {
  return new Promise((resolve, reject) => {
    const probe = net.createServer();
    probe.once('error', reject);
    probe.listen(0, '127.0.0.1', () => {
      const address = probe.address();
      if (typeof address !== 'object' || address === null) {
        probe.close();
        reject(new Error('Could not read a port from the probe socket'));
        return;
      }
      const { port } = address;
      probe.close(() => resolve(port));
    });
  });
};

// Only used when HERA_INTEG_PORT pins the port. Refuse to start if something is
// already there, rather than running the tests against a stranger's server.
const assertPortFree = (port: number): Promise<void> => {
  return new Promise((resolve, reject) => {
    const probe = net.createServer();
    probe.once('error', () => {
      reject(new Error(
        `HERA_INTEG_PORT=${port} is already in use. The integration tests write ` +
        `to the database they reach, so refusing to run against an existing server.`,
      ));
    });
    probe.listen(port, '127.0.0.1', () => { probe.close(() => resolve()); });
  });
};

const resolveServerPort = async (): Promise<number> => {
  const pinned = Number(process.env.HERA_INTEG_PORT);
  if (pinned) {
    await assertPortFree(pinned);
    return pinned;
  }
  return findFreePort();
};

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

export default async function setup() {
  const serverPort = await resolveServerPort();
  const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'hera-integ-'));
  const pyheraDir = path.join(tmpDir, '.pyhera');
  fs.mkdirSync(pyheraDir);

  fs.writeFileSync(path.join(pyheraDir, 'config.json'), JSON.stringify({
    [os.userInfo().username]: {
      dbIP: `127.0.0.1:${MONGO_PORT}`,
      dbName: 'olymp',
      username: 'hera',
      password: 'heracles',
    },
  }, null, 2));

  // Run the server from a scratch directory, not the repo. Test code creates
  // projects under os.getcwd(), so with cwd=PROJECT_ROOT every run littered the
  // working tree with a projects/ directory.
  const workDir = path.join(tmpDir, 'work');
  fs.mkdirSync(workDir);

  // server.py finds the client bundle from __file__ and its own directory goes on
  // sys.path, so neither depends on cwd.
  const proc = spawn('python', [
    path.join(PROJECT_ROOT, 'ui', 'server', 'server.py'),
    '--cors', 'all', '-y', '--jupyter-port', '0',
    '--port', String(serverPort),
  ], {
    cwd: workDir,
    env: { ...process.env, HOME: tmpDir },
    stdio: 'ignore',
  });

  const pid = proc.pid;
  if (!pid) {
    fs.rmSync(tmpDir, { recursive: true, force: true });
    throw new Error('Failed to spawn server process');
  }

  let exited = false;
  proc.on('exit', () => { exited = true; });

  let ready = false;
  const deadline = Date.now() + 30_000;
  while (Date.now() < deadline && !exited) {
    try {
      const r = await fetch(`http://localhost:${serverPort}/healthz`);
      if (r.ok) {
        ready = true;
        break;
      }
    } catch { /* not ready yet */ }
    await new Promise(r => setTimeout(r, 500));
  }

  if (!ready) {
    try { process.kill(pid, 'SIGTERM'); } catch { /* already exited */ }
    fs.rmSync(tmpDir, { recursive: true, force: true });
    throw new Error(
      exited
        ? `Server exited before it answered on port ${serverPort}`
        : `Server did not answer on port ${serverPort} within 30s`,
    );
  }

  console.log(`[globalSetup] Server ready on port ${serverPort} (mongo ${MONGO_PORT})`);
  fs.writeFileSync(PORT_FILE, String(serverPort));

  return function teardown() {
    try { process.kill(pid, 'SIGTERM'); } catch { /* already exited */ }
    try { fs.unlinkSync(PORT_FILE); } catch { /* ok */ }
  };
}
