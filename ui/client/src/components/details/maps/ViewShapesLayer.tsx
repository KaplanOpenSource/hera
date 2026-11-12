import { useEffect } from 'react';
import { useMap } from 'react-leaflet';
import L from 'leaflet';
import type { AnyProps, PointFeature } from 'supercluster';
import { serverApi } from '../../../api/api';
import { pointToDot } from './pointToDot';

const countInView = (points: PointFeature<AnyProps>[], view: L.LatLngBounds): number => {
  let count = 0;
  for (const point of points) {
    const [lng, lat] = point.geometry.coordinates;
    if (view.contains([lat, lng])) {
      count++;
    }
  }
  return count;
};

// Full shapes for the view, once it holds at most `limit` dots; the dots stay loaded.
export const ViewShapesLayer = ({
  path,
  points,
  limit,
  onShowing,
}: {
  path: string,
  points: PointFeature<AnyProps>[],
  limit: number,
  onShowing: (showing: boolean) => void,
}) => {
  const map = useMap();

  useEffect(() => {
    let shown: L.GeoJSON | null = null;
    // Each view gets a new ticket; a reply for an older ticket is dropped.
    let latest = 0;
    const clear = () => {
      if (shown) {
        map.removeLayer(shown);
        shown = null;
      }
    };
    const update = async () => {
      latest++;
      const ticket = latest;
      const view = map.getBounds();
      if (countInView(points, view) > limit) {
        clear();
        onShowing(false);
        return;
      }
      const answer = await serverApi.geoShapes(path, [view.getWest(), view.getSouth(), view.getEast(), view.getNorth()]);
      if (ticket !== latest || !answer) {
        return;
      }
      clear();
      shown = L.geoJSON(answer.shapes, { pointToLayer: pointToDot });
      map.addLayer(shown);
      onShowing(true);
    };
    update();
    map.on('moveend', update);
    return () => {
      latest++;
      map.off('moveend', update);
      clear();
    };
  }, [map, path, points, limit, onShowing]);

  return null;
};
