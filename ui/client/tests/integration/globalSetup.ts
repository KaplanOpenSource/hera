/// <reference types="node" />
import { spawn } from 'child_process';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';

const PORT_FILE = path.join(os.tmpdir(), 'hera-integ-port');
const SERVER_PORT = 8000;
const MONGO_PORT = 27018;

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

  const proc = spawn('python', [
    'ui/server/server.py', '--cors', 'all', '-y', '--jupyter-port', '0',
  ], {
    cwd: PROJECT_ROOT,
    env: { ...process.env, HOME: tmpDir },
    stdio: ['pipe', 'pipe', 'pipe'],
  });

  const pid = proc.pid;
  if (!pid) {
    fs.rmSync(tmpDir, { recursive: true, force: true });
    throw new Error('Failed to spawn server process');
  }

  const serverLog = path.join(tmpDir, 'server.log');
  const logStream = fs.createWriteStream(serverLog);
  proc.stdout?.pipe(logStream);
  proc.stderr?.pipe(logStream);

  let earlyExit: number | null = null;
  let stderr = '';
  proc.stderr?.on('data', (chunk: Buffer) => { stderr += chunk.toString(); });
  proc.on('exit', (code) => { earlyExit = code; });

  const deadline = Date.now() + 30_000;
  while (Date.now() < deadline) {
    if (earlyExit !== null) {
      fs.rmSync(tmpDir, { recursive: true, force: true });
      throw new Error(
        `Server exited early with code ${earlyExit}`
        + (stderr ? `\n${stderr.slice(-500)}` : ''),
      );
    }
    try {
      const r = await fetch(`http://localhost:${SERVER_PORT}/healthz`);
      if (r.ok) break;
    } catch { /* not ready yet */ }
    await new Promise(r => setTimeout(r, 500));
  }

  console.log(`[globalSetup] Server ready on port ${SERVER_PORT} (mongo ${MONGO_PORT})`);
  fs.writeFileSync(PORT_FILE, String(SERVER_PORT));

  return function teardown() {
    try { process.kill(pid, 'SIGTERM'); } catch { /* already exited */ }
    try { fs.unlinkSync(PORT_FILE); } catch { /* ok */ }
  };
}
