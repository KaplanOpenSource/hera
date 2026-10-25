import { describe, it, expect, vi, afterEach, beforeEach } from 'vitest';
import { createElement } from 'react';
import { render, waitFor } from '@testing-library/react';
import { useNodeCatalog, NodeCatalogReader } from '../src/components/workflow/useNodeCatalog';
import { mockApiCall } from './helpers/mockApiCall';

let getJson: ReturnType<typeof vi.fn>;

afterEach(() => {
  vi.restoreAllMocks();
});

// The store is a module singleton; reset it so each test starts empty.
beforeEach(() => {
  vi.clearAllMocks();
  ({ getJson } = mockApiCall());
  useNodeCatalog.setState({ catalog: [] });
});

describe('useNodeCatalog store', () => {
  it('loadCatalog fetches the catalog and stores it', async () => {
    getJson.mockResolvedValue([{ type: 'a.B', parameters: [] }]);

    await useNodeCatalog.getState().loadCatalog();

    expect(useNodeCatalog.getState().catalog).toHaveLength(1);
    expect(useNodeCatalog.getState().catalog[0].type).toBe('a.B');
    expect(getJson).toHaveBeenCalledTimes(1);
  });

  it('does not refetch once the catalog is loaded', async () => {
    getJson.mockResolvedValue([{ type: 'a.B', parameters: [] }]);

    await useNodeCatalog.getState().loadCatalog();
    await useNodeCatalog.getState().loadCatalog();

    expect(getJson).toHaveBeenCalledTimes(1);
  });

  it('leaves the catalog empty when the call fails', async () => {
    getJson.mockResolvedValue(null);

    await useNodeCatalog.getState().loadCatalog();

    expect(useNodeCatalog.getState().catalog).toEqual([]);
  });

  it('NodeCatalogReader triggers the fetch once on mount', async () => {
    getJson.mockResolvedValue([{ type: 'a.B', parameters: [] }]);

    render(createElement(NodeCatalogReader));

    await waitFor(() => expect(useNodeCatalog.getState().catalog).toHaveLength(1));
    expect(getJson).toHaveBeenCalledTimes(1);
  });
});
