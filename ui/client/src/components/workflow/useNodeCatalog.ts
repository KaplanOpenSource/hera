import { useEffect } from 'react';
import { create } from 'zustand';
import { BASEURL } from '../../shared/baseurl';
import { NodeCatalogEntry } from './nodeCatalog';

interface NodeCatalogStore {
  catalog: NodeCatalogEntry[];
  loadCatalog: () => Promise<void>;
}

// The Hermes node catalog (available node types + their parameters/outputs),
// fetched once from the server's /node-catalog endpoint and shared across the
// editor via this store — read it anywhere instead of fetching per node.
export const useNodeCatalog = create<NodeCatalogStore>((set, get) => ({
  catalog: [],
  loadCatalog: async () => {
    // Already loaded — the catalog doesn't change during a session.
    if (get().catalog.length > 0) {
      return;
    }
    const response = await fetch(`${BASEURL}/node-catalog`);
    if (!response.ok) {
      // The server's error middleware returns { error, traceback }; log the real
      // traceback instead of just the bare status so failures are debuggable here.
      const problem = await response.json() as { error?: string; traceback?: string };
      console.error(`node catalog failed (${response.status}): ${problem.error ?? response.statusText}\n${problem.traceback ?? ''}`);
      return;
    }
    set({ catalog: await response.json() as NodeCatalogEntry[] });
  },
}));

// Kicks off the one-time catalog fetch; render once inside the editor. A failed
// fetch must not crash the editor — log and leave the catalog empty.
export const NodeCatalogReader = () => {
  useEffect(() => {
    useNodeCatalog.getState().loadCatalog().catch(error => console.error('node catalog', error));
  }, []);
  return null;
};
