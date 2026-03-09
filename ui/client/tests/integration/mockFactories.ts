/** Shared mock factories for integration tests.
 *  Used via: vi.mock('...path...', async () => (await import('./mockFactories')).factoryName())
 */

export const createServerConstantsMock = () => {
  const state = { dataTypes: {}, readAllConstants: async () => {} };
  const useServerConstants = Object.assign(() => state, { getState: () => state });
  return { useServerConstants, ServerConstantReader: () => null };
};

export const createBaseurlMock = (port: number) => ({
  BASEURL: `http://localhost:${port}`,
});
