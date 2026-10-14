import { describe, it, expect, vi, beforeEach } from 'vitest';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFetch = vi.fn();
global.fetch = mockFetch;

const { fetchGeoJson } = await import('../src/io/fetchGeoJson');

const feature = { type: 'FeatureCollection', features: [] };

const jsonAnswer = (data: any) => ({
  ok: true,
  headers: { get: () => 'application/json' },
  json: () => Promise.resolve(data),
});

beforeEach(() => {
  vi.clearAllMocks();
});

describe('fetchGeoJson', () => {
  it('returns the GeoJSON the server read', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer(feature));

    expect(await fetchGeoJson('/data/a.geojson')).toEqual(feature);
  });

  it('sends the path encoded', async () => {
    mockFetch.mockResolvedValueOnce(jsonAnswer(feature));

    await fetchGeoJson('/my data/a b.geojson');

    expect(mockFetch.mock.calls[0][0]).toBe('http://test/geojson?path=%2Fmy%20data%2Fa%20b.geojson');
  });

  it('returns null when the server could not read it', async () => {
    mockFetch.mockResolvedValueOnce({
      ok: false,
      status: 404,
      headers: { get: () => 'application/json' },
      json: () => Promise.resolve({}),
    });

    expect(await fetchGeoJson('/missing.geojson')).toBeNull();
  });

  it('returns null when the server has no such route and answers HTML', async () => {
    mockFetch.mockResolvedValueOnce({
      ok: true,
      status: 200,
      headers: { get: () => 'text/html' },
      json: () => Promise.reject(new Error('not json')),
    });

    expect(await fetchGeoJson('/data/a.geojson')).toBeNull();
  });
});
