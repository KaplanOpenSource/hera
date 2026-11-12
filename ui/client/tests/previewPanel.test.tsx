import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, cleanup } from '@testing-library/react';

vi.mock('../src/shared/baseurl', () => ({ BASEURL: 'http://test' }));

// Every way out of the client, so a call by any route would show up here.
const mockFetch = vi.fn();
global.fetch = mockFetch;
import { mockApiCall } from './helpers/mockApiCall';

let mockGetJson: ReturnType<typeof vi.fn>;
let mockPostJson: ReturnType<typeof vi.fn>;

vi.mock('../src/components/details/ImagePreview', () => ({
  isImageResource: (r: unknown) => typeof r === 'string' && r.endsWith('.png'),
  ImagePreview: ({ path }: { path: string }) => <div data-testid="image">{path}</div>,
}));
vi.mock('../src/components/details/maps/GeoShapesPreview', () => ({
  isGeoShapesResource: () => false,
  GeoShapesPreview: () => null,
}));
vi.mock('../src/components/details/maps/GeoTiffPreview', () => ({
  isGeoTiffResource: (r: unknown) => typeof r === 'string' && r.endsWith('.tif'),
  isHgtResource: (r: unknown) => typeof r === 'string' && r.endsWith('.hgt'),
  TiffSource: { File: 'file', Hgt: 'hgt' },
  GeoTiffPreview: ({ path, source }: { path: string, source?: string }) => (
    <div data-testid="tiff" data-source={source ?? 'file'}>{path}</div>
  ),
}));
vi.mock('../src/components/details/maps/TileMapView', () => ({
  isTileUrl: () => false,
  TileMapView: () => null,
}));

const { PreviewPanel, hasPreview } = await import('../src/components/details/PreviewPanel');
const { useProjectStore } = await import('../src/stores/useProjectStore');

afterEach(() => { cleanup(); });

const doc = (oid: string, resource: string) => ({
  _id: { $oid: oid },
  _cls: 'Metadata.Cache',
  projectName: 'P',
  desc: { datasourceName: oid },
  type: 'T',
  resource,
  dataFormat: 'string',
});

beforeEach(() => {
  vi.clearAllMocks();
  const call = mockApiCall();
  mockGetJson = call.getJson;
  mockPostJson = call.postJson;
  useProjectStore.setState({
    currProjectName: 'P',
    currProject: { name: 'P', documents: [doc('doc1', '/pics/a.png')] as any },
    editedDocs: {},
  });
});

describe('PreviewPanel', () => {
  it('reads the resource from the store and makes no call', () => {
    const { getByTestId } = render(<PreviewPanel docid="doc1" />);

    expect(getByTestId('image').textContent).toBe('/pics/a.png');
    expect(mockFetch).not.toHaveBeenCalled();
    expect(mockGetJson).not.toHaveBeenCalled();
    expect(mockPostJson).not.toHaveBeenCalled();
  });

  it('shows nothing for a document the project does not have', () => {
    const { container } = render(<PreviewPanel docid="nope" />);

    expect(container.innerHTML).toBe('');
    expect(mockGetJson).not.toHaveBeenCalled();
  });
});

describe('PreviewPanel for a .hgt file', () => {
  it('has a preview', () => {
    expect(hasPreview(doc('hgt', '/data/N33E035.hgt') as any)).toBe(true);
  });

  it('shows it in the TIFF preview, read through the server', () => {
    useProjectStore.setState({
      currProject: { name: 'P', documents: [doc('hgt', '/data/N33E035.hgt')] as any },
    });

    const { getByTestId } = render(<PreviewPanel docid="hgt" />);

    expect(getByTestId('tiff').textContent).toBe('/data/N33E035.hgt');
    expect(getByTestId('tiff').dataset.source).toBe('hgt');
  });
});
