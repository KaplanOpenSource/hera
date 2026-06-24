import { describe, it, expect, vi, afterEach } from 'vitest';
import { renderHook, waitFor } from '@testing-library/react';
import { useNodeCatalog } from '../src/components/workflow/useNodeCatalog';

// The hook fetches the /node-catalog REST endpoint, so stub global fetch — never
// hit the real network.
afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe('useNodeCatalog', () => {
  it('fetches the catalog once and returns it', async () => {
    const fetchMock = vi.fn().mockResolvedValue({
      ok: true,
      json: async () => [{ type: 'a.B', parameters: [] }],
    });
    vi.stubGlobal('fetch', fetchMock);

    const { result } = renderHook(() => useNodeCatalog());
    await waitFor(() => expect(result.current.catalog.length).toBe(1));
    expect(result.current.catalog[0].type).toBe('a.B');
    expect(fetchMock).toHaveBeenCalledTimes(1);
  });

  it('leaves the catalog empty when the response is not ok', async () => {
    const fetchMock = vi.fn().mockResolvedValue({ ok: false, status: 500, statusText: 'err' });
    vi.spyOn(console, 'error').mockImplementation(() => {});
    vi.stubGlobal('fetch', fetchMock);

    const { result } = renderHook(() => useNodeCatalog());
    await waitFor(() => expect(fetchMock).toHaveBeenCalled());
    expect(result.current.catalog).toEqual([]);
  });
});
