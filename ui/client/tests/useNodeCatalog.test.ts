import { describe, it, expect, vi, afterEach, beforeEach } from 'vitest';
import { createElement } from 'react';
import { render, waitFor } from '@testing-library/react';
import { useNodeCatalog, NodeCatalogReader } from '../src/components/workflow/useNodeCatalog';

// The store fetches the /node-catalog REST endpoint, so stub global fetch — never
// hit the real network.
afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

// The store is a module singleton; reset it so each test starts empty.
beforeEach(() => {
  useNodeCatalog.setState({ catalog: [] });
});

describe('useNodeCatalog store', () => {
  it('loadCatalog fetches the catalog and stores it', async () => {
    const fetchMock = vi.fn().mockResolvedValue({
      ok: true,
      json: async () => [{ type: 'a.B', parameters: [] }],
    });
    vi.stubGlobal('fetch', fetchMock);

    await useNodeCatalog.getState().loadCatalog();

    expect(useNodeCatalog.getState().catalog).toHaveLength(1);
    expect(useNodeCatalog.getState().catalog[0].type).toBe('a.B');
    expect(fetchMock).toHaveBeenCalledTimes(1);
  });

  it('does not refetch once the catalog is loaded', async () => {
    const fetchMock = vi.fn().mockResolvedValue({
      ok: true,
      json: async () => [{ type: 'a.B', parameters: [] }],
    });
    vi.stubGlobal('fetch', fetchMock);

    await useNodeCatalog.getState().loadCatalog();
    await useNodeCatalog.getState().loadCatalog();

    expect(fetchMock).toHaveBeenCalledTimes(1);
  });

  it('leaves the catalog empty when the response is not ok', async () => {
    const fetchMock = vi.fn().mockResolvedValue({
      ok: false,
      status: 500,
      statusText: 'err',
      json: async () => ({ error: 'boom' }),
    });
    vi.spyOn(console, 'error').mockImplementation(() => {});
    vi.stubGlobal('fetch', fetchMock);

    await useNodeCatalog.getState().loadCatalog();

    expect(useNodeCatalog.getState().catalog).toEqual([]);
  });

  it('NodeCatalogReader triggers the fetch once on mount', async () => {
    const fetchMock = vi.fn().mockResolvedValue({
      ok: true,
      json: async () => [{ type: 'a.B', parameters: [] }],
    });
    vi.stubGlobal('fetch', fetchMock);

    render(createElement(NodeCatalogReader));

    await waitFor(() => expect(useNodeCatalog.getState().catalog).toHaveLength(1));
    expect(fetchMock).toHaveBeenCalledTimes(1);
  });
});
