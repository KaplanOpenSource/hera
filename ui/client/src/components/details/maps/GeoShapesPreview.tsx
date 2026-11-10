import { useEffect, useState } from 'react';
import { Box, CircularProgress, Typography, useTheme } from '@mui/material';
import { GeoJSON as GeoJSONLayer, MapContainer, TileLayer, useMap } from 'react-leaflet';
import L from 'leaflet';
import type { FeatureCollection } from 'geojson';
import type { AnyProps, PointFeature } from 'supercluster';
import 'leaflet/dist/leaflet.css';
import { serverApi } from '../../../api/api';
import { useProjectStore } from '../../../stores/useProjectStore';
import { ClusterLayer } from './ClusterLayer';
import { InvalidateOnResize } from './InvalidateOnResize';
import { pointToDot } from './pointToDot';
import { findConfigDoc } from './projectConfigDoc';
import { TileServerButton } from './TileServerButton';
import { darkMapControlsSx, TILE_ATTRIBUTION, TILE_URL } from './mapDarkStyles';

const GEO_SHAPES_EXTENSIONS = /\.(geojson|geo\.json|shp|fgb)$/i;

export const isGeoShapesResource = (resource: unknown): resource is string => {
  return typeof resource === 'string' && GEO_SHAPES_EXTENSIONS.test(resource);
};

// Above this many points, draw clusters instead of one layer per point.
export const CLUSTER_ABOVE = 1000;

// The points to cluster: only when every shape is a point and there are many.
export const pointsToCluster = (shapes: FeatureCollection): PointFeature<AnyProps>[] | null => {
  if (shapes.features.length <= CLUSTER_ABOVE) {
    return null;
  }
  for (const feature of shapes.features) {
    if (feature.geometry?.type !== 'Point') {
      return null;
    }
  }
  return shapes.features as PointFeature<AnyProps>[];
};

// Straight from the coordinates; L.geoJSON would build a layer per point.
const pointsBounds = (points: PointFeature<AnyProps>[]): L.LatLngBounds => {
  const bounds = L.latLngBounds([]);
  for (const point of points) {
    const [lng, lat] = point.geometry.coordinates;
    bounds.extend([lat, lng]);
  }
  return bounds;
};

type MapState = {
  shapes: FeatureCollection | null,
  points: PointFeature<AnyProps>[] | null,
  bounds: L.LatLngBounds | null,
  hasError: boolean,
};

const FitBounds = ({
  bounds,
}: {
  bounds: L.LatLngBounds,
}) => {
  const map = useMap();
  useEffect(() => {
    map.fitBounds(bounds, { padding: [20, 20] });
  }, [bounds, map]);
  return null;
};

export const GeoShapesPreview = ({
  path,
}: {
  path: string,
}) => {
  const [mapState, setMapState] = useState<MapState>({ shapes: null, points: null, bounds: null, hasError: false });
  const dark = useTheme().palette.mode === 'dark';
  const defaultTileServer = useProjectStore(state => findConfigDoc(state.currProject)?.desc?.defaultTileServer);
  let tileUrl = TILE_URL;
  let tileAttribution = TILE_ATTRIBUTION;
  if (defaultTileServer) {
    tileUrl = defaultTileServer;
    tileAttribution = '';
  }

  useEffect(() => {
    (async () => {
      setMapState({ shapes: null, points: null, bounds: null, hasError: false });
      const shapes = await serverApi.geoShapes(path);
      if (shapes) {
        const points = pointsToCluster(shapes);
        let bounds: L.LatLngBounds;
        if (points) {
          bounds = pointsBounds(points);
        } else {
          bounds = L.geoJSON(shapes).getBounds();
        }
        if (bounds.isValid()) {
          setMapState({ shapes, points, bounds, hasError: false });
          return;
        }
      }
      setMapState({ shapes: null, points: null, bounds: null, hasError: true });
    })();
  }, [path]);

  const centered = { display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' };

  return (
    mapState.hasError ? (
      <Box sx={centered}><Typography color="text.secondary">Map unavailable</Typography></Box>
    ) : !mapState.shapes || !mapState.bounds ? (
      <Box sx={centered}><CircularProgress /></Box>
    ) : (
      <Box sx={{ position: 'relative', height: '100%', width: '100%', ...(dark && darkMapControlsSx) }}>
        <Box sx={{ position: 'absolute', top: 10, right: 10, zIndex: 1000 }}>
          <TileServerButton currentUrl={tileUrl} />
        </Box>
        {/* Canvas, not SVG, so many shapes stay usable. */}
        <MapContainer
          center={[32, 35]}
          zoom={8}
          preferCanvas
          style={{ height: '100%', width: '100%' }}
        >
          <TileLayer url={tileUrl} attribution={tileAttribution} />
          {mapState.points && <ClusterLayer points={mapState.points} />}
          {!mapState.points && <GeoJSONLayer data={mapState.shapes} pointToLayer={pointToDot} />}
          <FitBounds bounds={mapState.bounds} />
          <InvalidateOnResize />
        </MapContainer>
      </Box>
    )
  );
};
