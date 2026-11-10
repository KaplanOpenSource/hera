import { useEffect, useMemo } from 'react';
import { useMap } from 'react-leaflet';
import L from 'leaflet';
import Supercluster, { type AnyProps, type PointFeature } from 'supercluster';
import { pointToDot } from './pointToDot';

// Near points within this many pixels merge into one cluster.
const CLUSTER_RADIUS = 60;
// Past this zoom every point shows on its own.
const MAX_CLUSTER_ZOOM = 16;
const MIN_CLUSTER_SIZE = 24;
const MAX_CLUSTER_SIZE = 48;

// Bigger clusters get bigger circles, on a log scale.
const clusterSize = (count: number): number => {
  const size = MIN_CLUSTER_SIZE + 6 * Math.log10(count);
  return Math.min(MAX_CLUSTER_SIZE, Math.round(size));
};

const clusterIcon = (count: number, label: string | number): L.DivIcon => {
  const size = clusterSize(count);
  const style = [
    `width:${size}px`, `height:${size}px`, `line-height:${size}px`,
    'border-radius:50%', 'text-align:center', 'font:12px sans-serif', 'color:#fff',
    'background:rgba(51,136,255,0.75)', 'border:2px solid rgba(255,255,255,0.8)', 'box-sizing:border-box',
  ].join(';');
  return L.divIcon({ html: `<div style="${style}">${label}</div>`, className: '', iconSize: [size, size] });
};

// Many points as clusters: only what is in view becomes a Leaflet layer.
export const ClusterLayer = ({
  points,
}: {
  points: PointFeature<AnyProps>[],
}) => {
  const map = useMap();
  const index = useMemo(() => {
    return new Supercluster({ radius: CLUSTER_RADIUS, maxZoom: MAX_CLUSTER_ZOOM }).load(points);
  }, [points]);

  useEffect(() => {
    const group = L.layerGroup();
    map.addLayer(group);
    const redraw = () => {
      group.clearLayers();
      const view = map.getBounds();
      const zoom = Math.round(map.getZoom());
      const items = index.getClusters([view.getWest(), view.getSouth(), view.getEast(), view.getNorth()], zoom);
      for (const item of items) {
        const [lng, lat] = item.geometry.coordinates;
        const latlng = L.latLng(lat, lng);
        if (!('cluster' in item.properties)) {
          group.addLayer(pointToDot(item, latlng));
          continue;
        }
        const { cluster_id, point_count, point_count_abbreviated } = item.properties;
        const marker = L.marker(latlng, { icon: clusterIcon(point_count, point_count_abbreviated) });
        marker.on('click', () => {
          map.setView(latlng, index.getClusterExpansionZoom(cluster_id));
        });
        group.addLayer(marker);
      }
    };
    redraw();
    map.on('moveend', redraw);
    return () => {
      map.off('moveend', redraw);
      map.removeLayer(group);
    };
  }, [map, index]);

  return null;
};
