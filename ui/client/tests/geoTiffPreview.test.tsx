import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, cleanup, screen } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFromUrl = vi.fn();
vi.mock('geotiff', () => ({ fromUrl: (url: string) => mockFromUrl(url) }));

const mockLayer = vi.fn();
vi.mock('../src/components/details/maps/GeoTiffGridLayer', () => ({
  GeoTiffGridLayer: class {
    constructor(_image: unknown, min: number, max: number, options: object) {
      mockLayer(min, max, options);
    }
    addTo() { return this; }
    remove() {}
  },
}));
vi.mock('react-leaflet', () => ({
  MapContainer: ({ children }: { children: any }) => <div data-testid="map">{children}</div>,
  useMap: () => ({ fitBounds: () => {} }),
}));
vi.mock('../src/components/details/maps/InvalidateOnResize', () => ({
  InvalidateOnResize: () => null,
}));

const { GeoTiffPreview, isGeoTiffResource } = await import('../src/components/details/maps/GeoTiffPreview');

const fakeImage = {
  getWidth: () => 2000,
  getHeight: () => 1000,
  getGDALNoData: () => 255,
  readRasters: vi.fn(async () => [new Uint8Array([2, 255, 7, 4])]),
};

afterEach(() => { cleanup(); });

beforeEach(() => {
  vi.clearAllMocks();
  vi.spyOn(console, 'error').mockImplementation(() => {});
  mockFromUrl.mockResolvedValue({ getImage: async () => fakeImage });
});

describe('isGeoTiffResource', () => {
  it('matches .tif and .tiff only', () => {
    expect(isGeoTiffResource('/a/b.tif')).toBe(true);
    expect(isGeoTiffResource('/a/b.TIFF')).toBe(true);
    expect(isGeoTiffResource('/a/b.png')).toBe(false);
    expect(isGeoTiffResource({})).toBe(false);
  });
});

describe('GeoTiffPreview', () => {
  it('reads the file through /file and draws it with its value range', async () => {
    render(<GeoTiffPreview path="/data/land.tif" />);

    await screen.findByTestId('map');
    expect(mockFromUrl).toHaveBeenCalledWith('http://test/file/data/land.tif');
    expect(fakeImage.readRasters).toHaveBeenCalledWith({ samples: [0], width: 512, height: 256 });
    expect(mockLayer).toHaveBeenCalledWith(2, 7, expect.objectContaining({ minZoom: -10 }));
  });

  it('shows a message when the file cannot be read', async () => {
    mockFromUrl.mockRejectedValue(new Error('404'));

    render(<GeoTiffPreview path="/missing.tif" />);

    expect(await screen.findByText('TIFF unavailable')).toBeTruthy();
  });
});
