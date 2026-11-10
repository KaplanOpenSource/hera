import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, cleanup, findByTestId, findByText } from '@testing-library/react';
import type { FeatureCollection } from 'geojson';
import { serialize } from 'flatgeobuf/geojson';
import { mockApiCall } from './helpers/mockApiCall';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

vi.mock('react-leaflet', () => ({
  MapContainer: ({ children }: { children: any }) => <div>{children}</div>,
  TileLayer: ({ url }: { url: string }) => <div data-testid="tiles">{url}</div>,
  GeoJSON: ({ data, pointToLayer }: { data: FeatureCollection, pointToLayer?: unknown }) => (
    <div data-testid="shapes" data-dots={String(Boolean(pointToLayer))}>{data.features.length}</div>
  ),
  useMap: () => ({ fitBounds: () => {}, invalidateSize: () => {} }),
}));
vi.mock('../src/components/details/maps/InvalidateOnResize', () => ({
  InvalidateOnResize: () => null,
}));
vi.mock('../src/components/details/maps/ClusterLayer', () => ({
  ClusterLayer: ({ points }: { points: unknown[] }) => <div data-testid="clusters">{points.length}</div>,
}));

const { GeoShapesPreview, CLUSTER_ABOVE, pointsToCluster } = await import('../src/components/details/maps/GeoShapesPreview');
const { TILE_URL } = await import('../src/components/details/maps/mapDarkStyles');
const { useProjectStore } = await import('../src/stores/useProjectStore');

const square: FeatureCollection = {
  type: 'FeatureCollection',
  features: [{
    type: 'Feature',
    properties: { name: 'a' },
    geometry: { type: 'Polygon', coordinates: [[[35, 32], [35.1, 32], [35.1, 32.1], [35, 32.1], [35, 32]]] },
  }],
};

const squareBytes = () => {
  const bytes = serialize(square, 4326);
  return bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
};

let mockGetBytes: ReturnType<typeof vi.fn>;

const configDoc = (desc: object) => ({
  _id: { $oid: 'cfg' },
  _cls: 'Metadata.Cache',
  projectName: 'P',
  desc,
  type: 'P__config__',
  resource: '',
  dataFormat: 'string',
});

const setDocs = (documents: object[]) => {
  useProjectStore.setState({
    currProjectName: 'P',
    currProject: { name: 'P', documents: documents as any },
    editedDocs: {},
  });
};

afterEach(() => { cleanup(); });

beforeEach(() => {
  vi.clearAllMocks();
  mockGetBytes = mockApiCall().getBytes;
  mockGetBytes.mockResolvedValue(squareBytes());
});

describe('GeoShapesPreview tiles', () => {
  it('uses the project defaultTileServer', async () => {
    const url = 'http://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}';
    setDocs([configDoc({ defaultTileServer: url })]);

    const { container } = render(<GeoShapesPreview path="/a.shp" />);

    expect((await findByTestId(container, 'tiles')).textContent).toBe(url);
  });

  it('falls back to OSM when the project has no defaultTileServer', async () => {
    setDocs([configDoc({})]);

    const { container } = render(<GeoShapesPreview path="/a.shp" />);

    expect((await findByTestId(container, 'tiles')).textContent).toBe(TILE_URL);
  });

  it('falls back to OSM when the project has no config document', async () => {
    setDocs([]);

    const { container } = render(<GeoShapesPreview path="/a.shp" />);

    expect((await findByTestId(container, 'tiles')).textContent).toBe(TILE_URL);
  });
});

describe('GeoShapesPreview shapes', () => {
  it('draws the decoded shapes on the map', async () => {
    setDocs([]);

    const { container } = render(<GeoShapesPreview path="/a.shp" />);

    const shapes = await findByTestId(container, 'shapes');
    expect(shapes.textContent).toBe('1');
    expect(shapes.dataset.dots).toBe('true');
    expect(mockGetBytes).toHaveBeenCalledWith('/geoshapes', { path: '/a.shp' }, { label: 'geoshapes' });
  });

  it('says the map is unavailable when the server could not read the file', async () => {
    setDocs([]);
    mockGetBytes.mockResolvedValue(null);

    const { container } = render(<GeoShapesPreview path="/a.shp" />);

    expect(await findByText(container, 'Map unavailable')).toBeTruthy();
  });
});

describe('isGeoShapesResource', () => {
  it('takes .shp, .fgb, .geojson and .geo.json', async () => {
    const { isGeoShapesResource } = await import('../src/components/details/maps/GeoShapesPreview');

    for (const path of ['/a.shp', '/a.SHP', '/a.fgb', '/a.geojson', '/a.geo.json']) {
      expect(isGeoShapesResource(path)).toBe(true);
    }
  });

  it('refuses other files and non-strings', async () => {
    const { isGeoShapesResource } = await import('../src/components/details/maps/GeoShapesPreview');

    for (const path of ['/a.gpkg', '/a.shx', '/a.dbf', '/a.tif', '/a.json', 5, null]) {
      expect(isGeoShapesResource(path)).toBe(false);
    }
  });
});

describe('pointToDot', () => {
  it('draws a point as a small circle, not a pin icon', async () => {
    const L = (await import('leaflet')).default;
    const { pointToDot } = await import('../src/components/details/maps/pointToDot');

    const layer = pointToDot(null, L.latLng(32, 35));

    expect(layer).toBeInstanceOf(L.CircleMarker);
    expect((layer as L.CircleMarker).getLatLng()).toEqual(L.latLng(32, 35));
  });
});

const pointsCollection = (count: number): FeatureCollection => {
  const features: FeatureCollection['features'] = [];
  for (let i = 0; i < count; i++) {
    features.push({ type: 'Feature', properties: {}, geometry: { type: 'Point', coordinates: [35 + i * 1e-5, 32] } });
  }
  return { type: 'FeatureCollection', features };
};

const asBuffer = (shapes: FeatureCollection) => {
  const bytes = serialize(shapes, 4326);
  return bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
};

describe('pointsToCluster', () => {
  it('clusters many points', () => {
    const shapes = pointsCollection(CLUSTER_ABOVE + 1);

    expect(pointsToCluster(shapes)).toHaveLength(CLUSTER_ABOVE + 1);
  });

  it('does not cluster a few points', () => {
    expect(pointsToCluster(pointsCollection(CLUSTER_ABOVE))).toBeNull();
  });

  it('does not cluster when any shape is not a point', () => {
    const shapes = pointsCollection(CLUSTER_ABOVE + 1);
    shapes.features.push(square.features[0]);

    expect(pointsToCluster(shapes)).toBeNull();
  });
});

describe('GeoShapesPreview clusters', () => {
  it('draws many points as clusters, not as one layer per point', async () => {
    setDocs([]);
    mockGetBytes.mockResolvedValue(asBuffer(pointsCollection(CLUSTER_ABOVE + 1)));

    const { container } = render(<GeoShapesPreview path="/big.shp" />);

    expect((await findByTestId(container, 'clusters')).textContent).toBe(String(CLUSTER_ABOVE + 1));
    expect(container.querySelector('[data-testid="shapes"]')).toBeNull();
  });

  it('draws a few points one by one', async () => {
    setDocs([]);
    mockGetBytes.mockResolvedValue(asBuffer(pointsCollection(10)));

    const { container } = render(<GeoShapesPreview path="/small.shp" />);

    expect((await findByTestId(container, 'shapes')).textContent).toBe('10');
    expect(container.querySelector('[data-testid="clusters"]')).toBeNull();
  });
});
