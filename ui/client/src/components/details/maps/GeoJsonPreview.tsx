import { useEffect, useState } from 'react';
import { Box, CircularProgress, Typography, useTheme } from '@mui/material';
import { GeoJSON as GeoJSONLayer, MapContainer, TileLayer, useMap } from 'react-leaflet';
import L from 'leaflet';
import type { GeoJSON } from 'geojson';
import 'leaflet/dist/leaflet.css';
import { serverApi } from '../../../api/api';
import { useProjectStore } from '../../../stores/useProjectStore';
import { InvalidateOnResize } from './InvalidateOnResize';
import { darkMapControlsSx, TILE_ATTRIBUTION, TILE_URL } from './mapDarkStyles';

const GEOJSON_EXTENSIONS = /\.(geojson|geo\.json)$/i;

export const isGeoJsonResource = (resource: unknown): resource is string => {
  return typeof resource === 'string' && GEOJSON_EXTENSIONS.test(resource);
};

type MapState = {
  geojson: GeoJSON | null,
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

export const GeoJsonPreview = ({
  path,
}: {
  path: string,
}) => {
  const [mapState, setMapState] = useState<MapState>({ geojson: null, bounds: null, hasError: false });
  const dark = useTheme().palette.mode === 'dark';
  const defaultTileServer = useProjectStore(state => {
    const project = state.currProject;
    const config = project?.documents.find(d => d.type === `${project.name}__config__`);
    return config?.desc?.defaultTileServer;
  });
  let tileUrl = TILE_URL;
  let tileAttribution = TILE_ATTRIBUTION;
  if (defaultTileServer) {
    tileUrl = defaultTileServer;
    tileAttribution = '';
  }

  useEffect(() => {
    (async () => {
      setMapState({ geojson: null, bounds: null, hasError: false });
      const geojson = await serverApi.geoJson(path);
      if (geojson) {
        const bounds = L.geoJSON(geojson as any).getBounds();
        if (bounds.isValid()) {
          setMapState({ geojson, bounds, hasError: false });
          return;
        }
      }
      setMapState({ geojson: null, bounds: null, hasError: true });
    })();
  }, [path]);

  const centered = { display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' };

  return (
    mapState.hasError ? (
      <Box sx={centered}><Typography color="text.secondary">Map unavailable</Typography></Box>
    ) : !mapState.geojson || !mapState.bounds ? (
      <Box sx={centered}><CircularProgress /></Box>
    ) : (
      <Box sx={{ height: '100%', width: '100%', ...(dark && darkMapControlsSx) }}>
        <MapContainer
          center={[32, 35]}
          zoom={8}
          style={{ height: '100%', width: '100%' }}
        >
          <TileLayer url={tileUrl} attribution={tileAttribution} />
          <GeoJSONLayer data={mapState.geojson as any} />
          <FitBounds bounds={mapState.bounds} />
          <InvalidateOnResize />
        </MapContainer>
      </Box>
    )
  );
};
