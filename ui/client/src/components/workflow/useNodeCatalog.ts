import { useEffect, useState } from 'react';
import { BASEURL } from '../../shared/baseurl';
import { NodeCatalogEntry } from './nodeCatalog';

// Fetches the Hermes node catalog (available node types + their parameters) once
// from the server's /node-catalog endpoint. Shared by every node in the editor —
// fetch high and pass down rather than per node.
export const useNodeCatalog = () => {
  const [catalog, setCatalog] = useState<NodeCatalogEntry[]>([]);

  useEffect(() => {
    const run = async () => {
      const response = await fetch(`${BASEURL}/node-catalog`);
      if (!response.ok) {
        console.error(`node catalog: ${response.status} ${response.statusText}`);
        return;
      }
      setCatalog(await response.json() as NodeCatalogEntry[]);
    };
    // A failed catalog fetch must not crash the editor — log and leave it empty.
    run().catch(error => console.error('node catalog', error));
  }, []);

  return { catalog };
};
