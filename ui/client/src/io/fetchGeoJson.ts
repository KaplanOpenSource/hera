import type { GeoJSON } from 'geojson';
import { BASEURL } from '../shared/baseurl';

// Reads a GeoJSON file from the server's disk, converted to WGS84.
// Returns null when the server could not read it.
export const fetchGeoJson = async (path: string): Promise<GeoJSON | null> => {
  // no-store: the answer must never come from the browser cache.
  const response = await fetch(`${BASEURL}/geojson?path=${encodeURIComponent(path)}`, { cache: 'no-store' });
  const contentType = response.headers.get('content-type');
  // A server without this route answers with the SPA's HTML, not an error.
  if (!response.ok || !contentType?.includes('application/json')) {
    console.error(`geojson: ${response.status} ${contentType} for ${path}`);
    return null;
  }
  return await response.json() as GeoJSON;
};
