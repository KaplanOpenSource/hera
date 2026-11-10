import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, cleanup, screen } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

const mockFromUrl = vi.fn();
const mockFromArrayBuffer = vi.fn();
vi.mock('geotiff', () => ({
  fromUrl: (url: string) => mockFromUrl(url),
  fromArrayBuffer: (bytes: ArrayBuffer) => mockFromArrayBuffer(bytes),
}));

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

const { GeoTiffPreview, isGeoTiffResource, isHgtResource, TiffSource } = await import('../src/components/details/maps/GeoTiffPreview');
const { serverApi } = await import('../src/api/api');

const fakeImage = {
  getWidth: () => 2000,
  getHeight: () => 1000,
  getGDALNoData: () => 255,
  readRasters: vi.fn(async () => [new Uint8Array([2, 255, 7, 4])]),
};

afterEach(() => { cleanup(); });

beforeEach(() => {
  vi.restoreAllMocks();
  vi.clearAllMocks();
  vi.spyOn(console, 'error').mockImplementation(() => {});
  mockFromUrl.mockResolvedValue({ getImage: async () => fakeImage });
  mockFromArrayBuffer.mockResolvedValue({ getImage: async () => fakeImage });
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

describe('isHgtResource', () => {
  it('matches .hgt only', () => {
    expect(isHgtResource('/a/N33E035.hgt')).toBe(true);
    expect(isHgtResource('/a/N33E035.HGT')).toBe(true);
    expect(isHgtResource('/a/b.tif')).toBe(false);
    expect(isHgtResource({})).toBe(false);
  });
});

describe('GeoTiffPreview of a .hgt file', () => {
  it('opens the GeoTIFF bytes the server made from it', async () => {
    const bytes = new ArrayBuffer(8);
    const hgt = vi.spyOn(serverApi, 'hgt').mockResolvedValue(bytes);

    render(<GeoTiffPreview path="/data/N33E035.hgt" source={TiffSource.Hgt} />);

    await screen.findByTestId('map');
    expect(hgt).toHaveBeenCalledWith('/data/N33E035.hgt');
    expect(mockFromArrayBuffer).toHaveBeenCalledWith(bytes);
    expect(mockFromUrl).not.toHaveBeenCalled();
    expect(mockLayer).toHaveBeenCalledWith(2, 7, expect.objectContaining({ minZoom: -10 }));
  });

  it('shows a message when the server could not read it', async () => {
    vi.spyOn(serverApi, 'hgt').mockResolvedValue(null);

    render(<GeoTiffPreview path="/missing.hgt" source={TiffSource.Hgt} />);

    expect(await screen.findByText('TIFF unavailable')).toBeTruthy();
  });
});
