import { readFileSync } from 'fs';
import { join } from 'path';
import { tmpdir } from 'os';

// The port globalSetup started the server on. It picks a free port at run time,
// so it must be read back here, never hardcoded on this side.
const readSharedPort = (): number => {
  const pinned = Number(process.env.HERA_INTEG_PORT);
  if (pinned) {
    return pinned;
  }
  const parsed = parseInt(readFileSync(join(tmpdir(), 'hera-integ-port'), 'utf-8').trim(), 10);
  if (isNaN(parsed)) {
    throw new Error('hera-integ-port does not hold a port number');
  }
  return parsed;
};

export const SHARED_PORT = readSharedPort();
export const SHARED_SERVER_URL = `http://localhost:${SHARED_PORT}`;

export const createServerConstantsMock = () => {
  const state = { dataTypes: { STRING: 'string', JSON: 'json', JSON_DICT: 'JSON_DICT' }, readAllConstants: async () => {} };
  const useServerConstants = Object.assign(() => state, { getState: () => state });
  return { useServerConstants, ServerConstantReader: () => null };
};

export const createBaseurlMock = (port: number = SHARED_PORT) => ({
  BASEURL: `http://localhost:${port}`,
});

let keyCounter = 0;
export const createSnackbarMock = () => ({
  pushRunning: (_label: string) => `mock-key-${++keyCounter}`,
  pushError: (_message: string) => `mock-key-${++keyCounter}`,
  dismiss: (_key: unknown) => {},
});
