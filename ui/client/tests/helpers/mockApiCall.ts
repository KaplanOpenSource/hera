import { vi } from 'vitest';
import { ApiCall } from '../../src/api/ApiCall';

// Stubs the methods every api class calls, so a test can watch the paths
// and bodies sent to the server without a real fetch. getJson, postJson, getBytes
// and postBytes are protected, which is a compile time word only, so the spy reaches them.
export const mockApiCall = () => {
  const getJson = vi.fn();
  const postJson = vi.fn();
  const getBytes = vi.fn();
  const postBytes = vi.fn();
  vi.spyOn(ApiCall.prototype as any, 'getJson').mockImplementation(getJson);
  vi.spyOn(ApiCall.prototype as any, 'postJson').mockImplementation(postJson);
  vi.spyOn(ApiCall.prototype as any, 'getBytes').mockImplementation(getBytes);
  vi.spyOn(ApiCall.prototype as any, 'postBytes').mockImplementation(postBytes);
  return { getJson, postJson, getBytes, postBytes };
};
