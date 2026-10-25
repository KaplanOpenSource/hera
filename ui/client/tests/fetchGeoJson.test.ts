import { describe, it, expect, vi, beforeEach } from 'vitest';
import { mockApiCall } from './helpers/mockApiCall';

const { serverApi } = await import('../src/api/api');

describe('serverApi.geoJson', () => {
  let getJson: ReturnType<typeof vi.fn>;

  beforeEach(() => {
    vi.clearAllMocks();
    ({ getJson } = mockApiCall());
  });

  it('returns the GeoJSON the server read', async () => {
    const feature = { type: 'FeatureCollection', features: [] };
    getJson.mockResolvedValueOnce(feature);

    expect(await serverApi.geoJson('/data/a.geojson')).toEqual(feature);
  });

  it('sends the path, unencoded, for the base class to encode', async () => {
    getJson.mockResolvedValueOnce({ type: 'FeatureCollection', features: [] });

    await serverApi.geoJson('/my data/a b.geojson');

    expect(getJson).toHaveBeenCalledWith('/geojson', { path: '/my data/a b.geojson' }, { label: 'geojson' });
  });

  it('returns null when the server could not read it', async () => {
    getJson.mockResolvedValueOnce(null);

    expect(await serverApi.geoJson('/missing.geojson')).toBeNull();
  });
});
