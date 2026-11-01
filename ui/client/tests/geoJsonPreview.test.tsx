import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, cleanup, findByTestId } from '@testing-library/react';
import { mockApiCall } from './helpers/mockApiCall';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

vi.mock('react-leaflet', () => ({
  MapContainer: ({ children }: { children: any }) => <div>{children}</div>,
  TileLayer: ({ url }: { url: string }) => <div data-testid="tiles">{url}</div>,
  GeoJSON: () => null,
  useMap: () => ({ fitBounds: () => {}, invalidateSize: () => {} }),
}));
vi.mock('../src/components/details/maps/InvalidateOnResize', () => ({
  InvalidateOnResize: () => null,
}));

const { GeoJsonPreview } = await import('../src/components/details/maps/GeoJsonPreview');
const { TILE_URL } = await import('../src/components/details/maps/mapDarkStyles');
const { useProjectStore } = await import('../src/stores/useProjectStore');

const square = {
  type: 'Polygon',
  coordinates: [[[35, 32], [35.1, 32], [35.1, 32.1], [35, 32.1], [35, 32]]],
};

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
  const { getJson } = mockApiCall();
  getJson.mockResolvedValue(square);
});

describe('GeoJsonPreview tiles', () => {
  it('uses the project defaultTileServer', async () => {
    const url = 'http://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}';
    setDocs([configDoc({ defaultTileServer: url })]);

    const { container } = render(<GeoJsonPreview path="/a.geojson" />);

    expect((await findByTestId(container, 'tiles')).textContent).toBe(url);
  });

  it('falls back to OSM when the project has no defaultTileServer', async () => {
    setDocs([configDoc({})]);

    const { container } = render(<GeoJsonPreview path="/a.geojson" />);

    expect((await findByTestId(container, 'tiles')).textContent).toBe(TILE_URL);
  });

  it('falls back to OSM when the project has no config document', async () => {
    setDocs([]);

    const { container } = render(<GeoJsonPreview path="/a.geojson" />);

    expect((await findByTestId(container, 'tiles')).textContent).toBe(TILE_URL);
  });
});
