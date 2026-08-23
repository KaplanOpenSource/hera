/// <reference types="node" />
import { spawn } from 'child_process';
import * as fs from 'fs';
import * as net from 'net';
import * as os from 'os';
import * as path from 'path';

const PORT_FILE = path.join(os.tmpdir(), 'hera-integ-port');
const MONGO_PORT = 27018;

// The tests read the port back from PORT_FILE (see mockFactories.ts), so any free
// port works. Asking the OS for one keeps the run from colliding with whatever else
// listens on this machine — server.py would otherwise silently move to the next free
// port and nothing would ever answer where the tests look.
const findFreePort = (): Promise<number> => new Promise((resolve, reject) => {
  const probe = net.createServer();
  probe.on('error', reject);
  probe.listen(0, '127.0.0.1', () => {
    const { port } = probe.address() as net.AddressInfo;
    probe.close(() => resolve(port));
  });
});

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
  const SERVER_PORT = await findFreePort();
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

  fs.writeFileSync(PORT_FILE, String(SERVER_PORT));

  const proc = spawn('python', [
    'ui/server/server.py', '--cors', 'all', '-y', '--jupyter-port', '0',
    '--port', String(SERVER_PORT),
  ], {
    cwd: PROJECT_ROOT,
    env: { ...process.env, HOME: tmpDir },
    stdio: 'ignore',
  });

  const pid = proc.pid;
  if (!pid) {
    fs.rmSync(tmpDir, { recursive: true, force: true });
    throw new Error('Failed to spawn server process');
  }

  // /healthz answers as soon as uvicorn binds, but /exec replies WARMING_UP until the
  // background thread finished importing hera. Wait for the import, or the first tests
  // race it. The import takes tens of seconds on a cold CI runner.
  const deadline = Date.now() + 180_000;
  let ready = false;
  while (!ready && Date.now() < deadline) {
    try {
      const r = await fetch(`http://localhost:${SERVER_PORT}/ready`);
      if (r.ok) ready = Boolean((await r.json()).ready);
    } catch { /* not up yet */ }
    if (!ready) await new Promise(r => setTimeout(r, 500));
  }
  if (!ready) {
    try { process.kill(pid, 'SIGTERM'); } catch { /* already exited */ }
    fs.rmSync(tmpDir, { recursive: true, force: true });
    throw new Error(`Server did not become ready on port ${SERVER_PORT} within 180s`);
  }

  console.log(`[globalSetup] Server ready on port ${SERVER_PORT} (mongo ${MONGO_PORT})`);

  return function teardown() {
    try { process.kill(pid, 'SIGTERM'); } catch { /* already exited */ }
    try { fs.unlinkSync(PORT_FILE); } catch { /* ok */ }
  };
}
