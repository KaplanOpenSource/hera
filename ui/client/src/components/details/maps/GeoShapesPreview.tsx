import { useEffect, useState } from 'react';
import { Box, CircularProgress, Typography, useTheme } from '@mui/material';
import { GeoJSON as GeoJSONLayer, MapContainer, TileLayer, useMap } from 'react-leaflet';
import L from 'leaflet';
import type { FeatureCollection } from 'geojson';
import 'leaflet/dist/leaflet.css';
import { serverApi } from '../../../api/api';
import { useProjectStore } from '../../../stores/useProjectStore';
import { InvalidateOnResize } from './InvalidateOnResize';
import { findConfigDoc } from './projectConfigDoc';
import { TileServerButton } from './TileServerButton';
import { darkMapControlsSx, TILE_ATTRIBUTION, TILE_URL } from './mapDarkStyles';

const GEO_SHAPES_EXTENSIONS = /\.(geojson|geo\.json|shp|fgb)$/i;

export const isGeoShapesResource = (resource: unknown): resource is string => {
  return typeof resource === 'string' && GEO_SHAPES_EXTENSIONS.test(resource);
};

const DOT_RADIUS = 3;

// Small canvas circles, not Leaflet's pin icons, so many points stay usable.
export const pointToDot = (_feature: unknown, latlng: L.LatLng): L.Layer => {
  return L.circleMarker(latlng, { radius: DOT_RADIUS });
};

type MapState = {
  shapes: FeatureCollection | null,
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
  const [mapState, setMapState] = useState<MapState>({ shapes: null, bounds: null, hasError: false });
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
      setMapState({ shapes: null, bounds: null, hasError: false });
      const shapes = await serverApi.geoShapes(path);
      if (shapes) {
        const bounds = L.geoJSON(shapes).getBounds();
        if (bounds.isValid()) {
          setMapState({ shapes, bounds, hasError: false });
          return;
        }
      }
      setMapState({ shapes: null, bounds: null, hasError: true });
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
          <GeoJSONLayer data={mapState.shapes} pointToLayer={pointToDot} />
          <FitBounds bounds={mapState.bounds} />
          <InvalidateOnResize />
        </MapContainer>
      </Box>
    )
  );
};
