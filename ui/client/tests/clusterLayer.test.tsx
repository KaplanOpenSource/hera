import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, cleanup } from '@testing-library/react';
import L from 'leaflet';
import type { AnyProps, PointFeature } from 'supercluster';

// A map with a fixed view; layers added to it are only recorded.
const fakeMap = {
  view: L.latLngBounds([31, 34], [33, 36]),
  zoom: 5,
  layers: [] as L.Layer[],
  handlers: [] as (() => void)[],
  getBounds() { return this.view; },
  getZoom() { return this.zoom; },
  addLayer(layer: L.Layer) { this.layers.push(layer); },
  removeLayer(layer: L.Layer) { this.layers = this.layers.filter(l => l !== layer); },
  on(_name: string, handler: () => void) { this.handlers.push(handler); },
  off(_name: string, handler: () => void) { this.handlers = this.handlers.filter(h => h !== handler); },
  setView: vi.fn(),
};

vi.mock('react-leaflet', () => ({ useMap: () => fakeMap }));

const { ClusterLayer } = await import('../src/components/details/maps/ClusterLayer');

const pointsAt = (count: number, lng: number, lat: number): PointFeature<AnyProps>[] => {
  const points: PointFeature<AnyProps>[] = [];
  for (let i = 0; i < count; i++) {
    points.push({ type: 'Feature', properties: {}, geometry: { type: 'Point', coordinates: [lng + i * 1e-6, lat] } });
  }
  return points;
};

const drawn = (): L.Layer[] => {
  return (fakeMap.layers[0] as L.LayerGroup).getLayers();
};

afterEach(() => { cleanup(); });

beforeEach(() => {
  vi.clearAllMocks();
  fakeMap.zoom = 5;
  fakeMap.layers = [];
  fakeMap.handlers = [];
});

describe('ClusterLayer', () => {
  it('merges near points into one cluster with their count', () => {
    render(<ClusterLayer points={pointsAt(1500, 35, 32)} />);

    const layers = drawn();
    expect(layers).toHaveLength(1);
    expect(layers[0]).toBeInstanceOf(L.Marker);
    expect((layers[0] as L.Marker).getIcon().options.html).toContain('1.5k');
  });

  it('shows far points apart, as dots', () => {
    render(<ClusterLayer points={[...pointsAt(1, 34.5, 31.5), ...pointsAt(1, 35.5, 32.5)]} />);

    const layers = drawn();
    expect(layers).toHaveLength(2);
    expect(layers[0]).toBeInstanceOf(L.CircleMarker);
    expect(layers[1]).toBeInstanceOf(L.CircleMarker);
  });

  it('leaves out points outside the view', () => {
    render(<ClusterLayer points={[...pointsAt(1, 35, 32), ...pointsAt(1, 10, 10)]} />);

    expect(drawn()).toHaveLength(1);
  });

  it('splits clusters into dots when zoomed in and the map moves', () => {
    render(<ClusterLayer points={[...pointsAt(1, 35, 32), ...pointsAt(1, 35.001, 32)]} />);
    expect(drawn()).toHaveLength(1);

    fakeMap.zoom = 18;
    for (const handler of fakeMap.handlers) {
      handler();
    }

    expect(drawn()).toHaveLength(2);
  });

  it('zooms in on a cluster when it is clicked', () => {
    render(<ClusterLayer points={pointsAt(1500, 35, 32)} />);

    (drawn()[0] as L.Marker).fire('click');

    expect(fakeMap.setView).toHaveBeenCalledTimes(1);
    expect(fakeMap.setView.mock.calls[0][1]).toBeGreaterThan(fakeMap.zoom);
  });

  it('removes its layers and stops listening when it goes away', () => {
    const { unmount } = render(<ClusterLayer points={pointsAt(10, 35, 32)} />);

    unmount();

    expect(fakeMap.layers).toHaveLength(0);
    expect(fakeMap.handlers).toHaveLength(0);
  });

  it('draws nothing while hidden', () => {
    render(<ClusterLayer points={pointsAt(10, 35, 32)} hidden />);

    expect(fakeMap.layers).toHaveLength(0);
    expect(fakeMap.handlers).toHaveLength(0);
  });
});
