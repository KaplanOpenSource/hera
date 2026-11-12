import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, cleanup, waitFor } from '@testing-library/react';
import L from 'leaflet';
import type { FeatureCollection } from 'geojson';
import type { AnyProps, PointFeature } from 'supercluster';

// A map with a set view; layers added to it are only recorded.
const fakeMap = {
  view: L.latLngBounds([31.9, 34.9], [32.1, 35.1]),
  layers: [] as L.Layer[],
  handlers: [] as (() => void)[],
  getBounds() { return this.view; },
  addLayer(layer: L.Layer) { this.layers.push(layer); },
  removeLayer(layer: L.Layer) { this.layers = this.layers.filter(l => l !== layer); },
  on(_name: string, handler: () => void) { this.handlers.push(handler); },
  off(_name: string, handler: () => void) { this.handlers = this.handlers.filter(h => h !== handler); },
};

vi.mock('react-leaflet', () => ({ useMap: () => fakeMap }));

const { serverApi } = await import('../src/api/api');
const { ViewShapesLayer } = await import('../src/components/details/maps/ViewShapesLayer');

const pointsAt = (count: number, lng: number, lat: number): PointFeature<AnyProps>[] => {
  const points: PointFeature<AnyProps>[] = [];
  for (let i = 0; i < count; i++) {
    points.push({ type: 'Feature', properties: {}, geometry: { type: 'Point', coordinates: [lng + i * 1e-6, lat] } });
  }
  return points;
};

const squareNamed = (name: string): FeatureCollection => {
  return {
    type: 'FeatureCollection',
    features: [{
      type: 'Feature',
      properties: { name },
      geometry: { type: 'Polygon', coordinates: [[[35, 32], [35.01, 32], [35.01, 32.01], [35, 32.01], [35, 32]]] },
    }],
  };
};

const shownNames = (): string[] => {
  const names: string[] = [];
  for (const layer of fakeMap.layers) {
    (layer as L.GeoJSON).eachLayer(child => {
      names.push((child as L.Polygon).feature?.properties.name);
    });
  }
  return names;
};

const moveMap = () => {
  for (const handler of fakeMap.handlers) {
    handler();
  }
};

// Three dots in view, one far outside it.
const points = [...pointsAt(3, 35, 32), ...pointsAt(1, 10, 10)];

let geoShapes: ReturnType<typeof vi.spyOn>;
let onShowing: ReturnType<typeof vi.fn>;

afterEach(() => { cleanup(); });

beforeEach(() => {
  vi.restoreAllMocks();
  fakeMap.view = L.latLngBounds([31.9, 34.9], [32.1, 35.1]);
  fakeMap.layers = [];
  fakeMap.handlers = [];
  onShowing = vi.fn();
  geoShapes = vi.spyOn(serverApi, 'geoShapes').mockResolvedValue({ shapes: squareNamed('a'), dotsLimit: null });
});

describe('ViewShapesLayer', () => {
  it('asks for the full shapes of the view when it holds at most the limit', async () => {
    render(<ViewShapesLayer path="/big.shp" points={points} limit={3} onShowing={onShowing} />);

    await waitFor(() => { expect(onShowing).toHaveBeenCalledWith(true); });
    expect(geoShapes).toHaveBeenCalledWith('/big.shp', [34.9, 31.9, 35.1, 32.1]);
    expect(shownNames()).toEqual(['a']);
  });

  it('asks nothing and shows no shapes over the limit', async () => {
    render(<ViewShapesLayer path="/big.shp" points={points} limit={2} onShowing={onShowing} />);

    expect(onShowing).toHaveBeenCalledWith(false);
    expect(geoShapes).not.toHaveBeenCalled();
    expect(fakeMap.layers).toHaveLength(0);
  });

  it('drops the shapes when the view goes over the limit again', async () => {
    render(<ViewShapesLayer path="/big.shp" points={points} limit={3} onShowing={onShowing} />);
    await waitFor(() => { expect(shownNames()).toEqual(['a']); });

    fakeMap.view = L.latLngBounds([0, 0], [40, 40]);
    moveMap();

    expect(fakeMap.layers).toHaveLength(0);
    expect(onShowing).toHaveBeenLastCalledWith(false);
  });

  it('drops a reply for a view the map already left', async () => {
    let answerOld: (value: any) => void = () => {};
    geoShapes.mockImplementationOnce(() => new Promise(resolve => { answerOld = resolve; }));
    geoShapes.mockResolvedValueOnce({ shapes: squareNamed('new'), dotsLimit: null });
    render(<ViewShapesLayer path="/big.shp" points={points} limit={3} onShowing={onShowing} />);

    moveMap();
    await waitFor(() => { expect(shownNames()).toEqual(['new']); });
    answerOld({ shapes: squareNamed('old'), dotsLimit: null });
    await Promise.resolve();

    expect(shownNames()).toEqual(['new']);
  });

  it('keeps what it shows when a bbox call fails', async () => {
    render(<ViewShapesLayer path="/big.shp" points={points} limit={3} onShowing={onShowing} />);
    await waitFor(() => { expect(shownNames()).toEqual(['a']); });

    geoShapes.mockResolvedValueOnce(null);
    moveMap();
    await Promise.resolve();

    expect(shownNames()).toEqual(['a']);
  });

  it('removes its shapes and stops listening when it goes away', async () => {
    const { unmount } = render(<ViewShapesLayer path="/big.shp" points={points} limit={3} onShowing={onShowing} />);
    await waitFor(() => { expect(shownNames()).toEqual(['a']); });

    unmount();

    expect(fakeMap.layers).toHaveLength(0);
    expect(fakeMap.handlers).toHaveLength(0);
  });
});
