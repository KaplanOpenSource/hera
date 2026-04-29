import { startDockerEnv, type DockerEnv } from './dockerSetup';

export default async function setup() {
  const env: DockerEnv = await startDockerEnv({
    network: 'hera-test-shared-net',
    mongoContainer: 'hera-test-shared-mongo',
    serverContainer: 'hera-test-shared-server',
    serverPort: 8001,
    dbName: 'hera_test',
  });

  return function teardown() {
    env.cleanup();
  };
}
