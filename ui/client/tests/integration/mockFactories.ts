/** Shared mock factories for integration tests.
 *  Used via: vi.mock('...path...', async () => (await import('./mockFactories')).factoryName())
 */

export const createServerConstantsMock = () => {
  const state = { dataTypes: { STRING: 'string', JSON: 'json', JSON_DICT: 'JSON_DICT' }, readAllConstants: async () => {} };
  const useServerConstants = Object.assign(() => state, { getState: () => state });
  return { useServerConstants, ServerConstantReader: () => null };
};

export const createBaseurlMock = (port: number) => ({
  BASEURL: `http://localhost:${port}`,
});

let keyCounter = 0;
export const createSnackbarMock = () => ({
  pushRunning: (_label: string) => `mock-key-${++keyCounter}`,
  pushError: (_message: string) => `mock-key-${++keyCounter}`,
  dismiss: (_key: unknown) => {},
});
