/// <reference types="node" />
import { execSync } from 'child_process';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';

const PROJECT_ROOT = path.resolve(__dirname, '../../../../');

const docker = (cmd: string) =>
  execSync(`docker ${cmd}`, { encoding: 'utf-8', stdio: ['pipe', 'pipe', 'pipe'] }).trim();

export type DockerEnv = {
  mongoContainerId: string;
  serverContainerId: string;
  network: string;
  port: number;
  tmpDir: string;
};

export const checkDocker = (): boolean => {
  try {
    docker('info');
    docker('image inspect hera-server');
    return true;
  } catch {
    return false;
  }
};

export const startDockerEnv = async ({
  network,
  mongoContainer,
  serverContainer,
  serverPort,
  dbName,
}: {
  network: string;
  mongoContainer: string;
  serverContainer: string;
  serverPort: number;
  dbName: string;
}): Promise<DockerEnv> => {
  docker(`network create ${network}`);

  const mongoInitDir = path.join(PROJECT_ROOT, 'mongo-init.d');
  const mongoContainerId = docker([
    `run -d --name ${mongoContainer}`,
    `--network ${network}`,
    `-v ${mongoInitDir}:/docker-entrypoint-initdb.d`,
    `mongo:5.0`,
  ].join(' '));

  // Wait for MongoDB to be ready
  const mongoDeadline = Date.now() + 15_000;
  let mongoReady = false;
  while (Date.now() < mongoDeadline) {
    try {
      const out = docker(`exec ${mongoContainerId} mongosh --eval "db.runCommand({ping:1})" --quiet`);
      if (out.includes('ok')) { mongoReady = true; break; }
    } catch { /* not ready */ }
    await new Promise(r => setTimeout(r, 300));
  }
  if (!mongoReady) {
    stopDockerEnv({ mongoContainerId, serverContainerId: '', network, port: serverPort, tmpDir: '' });
    throw new Error(`MongoDB container "${mongoContainer}" did not become ready within 15 s`);
  }

  const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'hera-integ-'));
  const pyheraDir = path.join(tmpDir, '.pyhera');
  fs.mkdirSync(pyheraDir);

  fs.writeFileSync(path.join(pyheraDir, 'config.json'), JSON.stringify({
    root: {
      dbName,
      dbIP: mongoContainer,
      username: 'hera',
      password: 'heracles',
    },
  }, null, 2));

  const serverContainerId = docker([
    `run -d --name ${serverContainer}`,
    `--network ${network}`,
    `-p 127.0.0.1:${serverPort}:8000`,
    `-v ${PROJECT_ROOT}:/app`,
    `-v ${pyheraDir}:/root/.pyhera`,
    `hera-server`,
    `python ui/server/server.py --cors all -y --jupyter-port 0`,
  ].join(' '));

  // Wait for server to be ready
  const serverDeadline = Date.now() + 20_000;
  let serverReady = false;
  while (Date.now() < serverDeadline) {
    try {
      const r = await fetch(`http://localhost:${serverPort}/healthz`);
      if (r.ok) { serverReady = true; break; }
    } catch { /* not ready yet */ }
    await new Promise(r => setTimeout(r, 300));
  }
  if (!serverReady) {
    let logs = '';
    try { logs = docker(`logs --tail 30 ${serverContainerId} 2>&1`); } catch { /* can't get logs */ }
    stopDockerEnv({ mongoContainerId, serverContainerId, network, port: serverPort, tmpDir });
    throw new Error(
      `Server container "${serverContainer}" did not become ready within 20 s.\n`
      + (logs ? `Last 30 lines of container logs:\n${logs}` : 'Could not retrieve container logs.'),
    );
  }

  return { mongoContainerId, serverContainerId, network, port: serverPort, tmpDir };
};

export const stopDockerEnv = (env: DockerEnv) => {
  if (env.serverContainerId) {
    try { docker(`rm -f ${env.serverContainerId}`); } catch { /* already removed */ }
  }
  if (env.mongoContainerId) {
    try { docker(`rm -f ${env.mongoContainerId}`); } catch { /* already removed */ }
  }
  if (env.network) {
    try { docker(`network rm ${env.network}`); } catch { /* already removed */ }
  }
  if (env.tmpDir) {
    try { fs.rmSync(env.tmpDir, { recursive: true, force: true }); } catch { /* best effort */ }
  }
};
