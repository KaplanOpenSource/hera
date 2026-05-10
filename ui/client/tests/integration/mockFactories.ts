import { readFileSync } from 'fs';
import { join } from 'path';
import { tmpdir } from 'os';

let detectedPort = 8001;
try {
  const parsed = parseInt(readFileSync(join(tmpdir(), 'hera-integ-port'), 'utf-8').trim(), 10);
  if (!isNaN(parsed)) detectedPort = parsed;
} catch { /* port file not written yet — fall back to default */ }

export const SHARED_PORT = detectedPort;
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
