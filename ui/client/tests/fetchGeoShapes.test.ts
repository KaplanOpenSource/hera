import { describe, it, expect, vi, beforeEach } from 'vitest';
import type { FeatureCollection } from 'geojson';
import { readFileSync } from 'fs';
import { dirname, join } from 'path';
import { fileURLToPath } from 'url';
import { serialize } from 'flatgeobuf/geojson';
import { mockApiCall } from './helpers/mockApiCall';

const { serverApi } = await import('../src/api/api');

const shapes: FeatureCollection = {
  type: 'FeatureCollection',
  features: [
    {
      type: 'Feature',
      properties: { name: 'a', height: 10.5 },
      geometry: { type: 'Polygon', coordinates: [[[35, 32], [35.1, 32], [35.1, 32.1], [35, 32.1], [35, 32]]] },
    },
    {
      type: 'Feature',
      properties: { name: 'b', height: 21 },
      geometry: { type: 'Polygon', coordinates: [[[35.2, 32], [35.3, 32], [35.3, 32.1], [35.2, 32.1], [35.2, 32]]] },
    },
  ],
};

const DOTS_FIXTURE = join(dirname(fileURLToPath(import.meta.url)), 'fixtures', 'dots.fgb');

const asBuffer = (bytes: Uint8Array): ArrayBuffer => {
  return bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) as ArrayBuffer;
};

describe('serverApi.geoShapes', () => {
  let getBytes: ReturnType<typeof vi.fn>;

  beforeEach(() => {
    vi.clearAllMocks();
    ({ getBytes } = mockApiCall());
  });

  it('decodes the FlatGeobuf the server sent into the same shapes', async () => {
    getBytes.mockResolvedValueOnce(asBuffer(serialize(shapes, 4326)));

    const result = await serverApi.geoShapes('/data/a.shp');

    expect(result?.shapes.features.map(f => f.geometry)).toEqual(shapes.features.map(f => f.geometry));
    expect(result?.shapes.features.map(f => f.properties)).toEqual(shapes.features.map(f => f.properties));
    expect(result?.dotsLimit).toBeNull();
  });

  it('sends the path, unencoded, for the base class to encode', async () => {
    getBytes.mockResolvedValueOnce(asBuffer(serialize(shapes, 4326)));

    await serverApi.geoShapes('/my data/a b.shp');

    expect(getBytes).toHaveBeenCalledWith('/geoshapes', { path: '/my data/a b.shp' }, { label: 'geoshapes', silent: false });
  });

  it('returns null when the server could not read it', async () => {
    getBytes.mockResolvedValueOnce(null);

    expect(await serverApi.geoShapes('/missing.shp')).toBeNull();
  });

  it('reads the dots flag and its limit from the title', async () => {
    getBytes.mockResolvedValueOnce(asBuffer(new Uint8Array(readFileSync(DOTS_FIXTURE))));

    const result = await serverApi.geoShapes('/data/big.shp');

    expect(result?.dotsLimit).toBe(2);
    expect(result?.shapes.features).toHaveLength(3);
  });

  it('sends a bbox quietly, as west,south,east,north', async () => {
    getBytes.mockResolvedValueOnce(asBuffer(serialize(shapes, 4326)));

    await serverApi.geoShapes('/data/a.shp', [34.5, 31.5, 35.5, 32.5]);

    expect(getBytes).toHaveBeenCalledWith(
      '/geoshapes', { path: '/data/a.shp', bbox: '34.5,31.5,35.5,32.5' }, { label: 'geoshapes', silent: true },
    );
  });

  it('reads no bytes as no shapes', async () => {
    getBytes.mockResolvedValueOnce(new ArrayBuffer(0));

    expect(await serverApi.geoShapes('/data/a.shp', [0, 0, 1, 1])).toEqual({
      shapes: { type: 'FeatureCollection', features: [] },
      dotsLimit: null,
    });
  });
});
