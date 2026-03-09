/// <reference types="node" />
import { execSync } from 'child_process';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';

const PROJECT_ROOT = path.resolve(__dirname, '../../../../');

const assertDockerImagesExist = () => {
  try {
    execSync('docker image inspect hera-server', { stdio: ['pipe', 'pipe', 'pipe'] });
  } catch {
    throw new Error(
      'Required Docker image "hera-server" not found.\n'
      + 'Build it from the project root with:\n'
      + '  sh hera/scripts/docker_build.sh',
    );
  }
};

export type DockerEnv = {
  serverUrl: string;
  tmpDir: string;
  cleanup: () => void;
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
  assertDockerImagesExist();

  const docker = (cmd: string) => execSync(`docker ${cmd}`, { encoding: 'utf-8' }).trim();

  const cleanupDocker = () => {
    for (const name of [serverContainer, mongoContainer]) {
      try { docker(`rm -f ${name}`); } catch { /* ignore */ }
    }
    try { docker(`network rm ${network}`); } catch { /* ignore */ }
  };

  // Clean up leftovers from a previous failed run
  cleanupDocker();

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

  docker(`network create ${network}`);

  const mongoInitDir = path.join(PROJECT_ROOT, 'mongo-init.d');
  docker([
    `run -d --name ${mongoContainer}`,
    `--network ${network}`,
    `-v ${mongoInitDir}:/docker-entrypoint-initdb.d`,
    `mongo:5.0`,
  ].join(' '));

  // Wait for MongoDB to be ready
  const mongoStart = Date.now();
  while (Date.now() - mongoStart < 15000) {
    try {
      const out = docker(`exec ${mongoContainer} mongosh --eval "db.runCommand({ping:1})" --quiet`);
      if (out.includes('ok')) break;
    } catch { /* not ready */ }
    await new Promise(r => setTimeout(r, 300));
  }

  docker([
    `run -d --name ${serverContainer}`,
    `--network ${network}`,
    `-p 127.0.0.1:${serverPort}:8000`,
    `-v ${PROJECT_ROOT}:/app`,
    `-v ${pyheraDir}:/root/.pyhera`,
    `hera-server`,
    `python ui/server/server.py --cors`,
  ].join(' '));

  const serverUrl = `http://localhost:${serverPort}`;

  // Wait for server to be ready
  const start = Date.now();
  while (Date.now() - start < 20000) {
    try {
      const r = await fetch(`${serverUrl}/healthz`);
      if (r.ok) break;
    } catch { /* not ready yet */ }
    await new Promise(r => setTimeout(r, 300));
  }

  return {
    serverUrl,
    tmpDir,
    cleanup: () => {
      cleanupDocker();
      try { fs.rmSync(tmpDir, { recursive: true, force: true }); } catch { /* ignore */ }
    },
  };
};
