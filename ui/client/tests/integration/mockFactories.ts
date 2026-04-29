export const SHARED_PORT = 8001;
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
