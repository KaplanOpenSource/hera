import { useEffect, useState } from 'react';
import { Box, CircularProgress, Typography } from '@mui/material';
import { GeoJSON as GeoJSONLayer, MapContainer, TileLayer, useMap } from 'react-leaflet';
import L from 'leaflet';
import type { GeoJSON } from 'geojson';
import 'leaflet/dist/leaflet.css';
import { fetchPython } from '../../../io/fetchPython';
import { InvalidateOnResize } from './InvalidateOnResize';

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

const loadGeoJson = async (path: string): Promise<GeoJSON | null> => {
  const { data } = await fetchPython({
    results: ['geojson_data'],
    code: `
import geopandas as gpd
import json
gdf = gpd.read_file('${path}')
if gdf.crs and not gdf.crs.is_geographic:
    gdf = gdf.to_crs(epsg=4326)
else:
    bounds = gdf.total_bounds
    if abs(bounds[0]) > 360 or abs(bounds[1]) > 360:
        gdf = gdf.set_crs(epsg=2039, allow_override=True)
        gdf = gdf.to_crs(epsg=4326)
geojson_data = json.loads(gdf.to_json())
`,
  });
  return data?.geojson_data ?? null;
};

export const GeoJsonPreview = ({
  path,
}: {
  path: string,
}) => {
  const [state, setState] = useState<MapState>({ geojson: null, bounds: null, hasError: false });

  useEffect(() => {
    setState({ geojson: null, bounds: null, hasError: false });
    (async () => {
      const geojson = await loadGeoJson(path);
      if (!geojson) {
        setState({ geojson: null, bounds: null, hasError: true });
        return;
      }
      const bounds = L.geoJSON(geojson as any).getBounds();
      if (!bounds.isValid()) {
        setState({ geojson: null, bounds: null, hasError: true });
        return;
      }
      setState({ geojson, bounds, hasError: false });
    })();
  }, [path]);

  const centered = { display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' };

  return (
    state.hasError ? (
      <Box sx={centered}><Typography color="text.secondary">Map unavailable</Typography></Box>
    ) : !state.geojson || !state.bounds ? (
      <Box sx={centered}><CircularProgress /></Box>
    ) : (
      <MapContainer
        center={[32, 35]}
        zoom={8}
        style={{ height: '100%', width: '100%' }}
      >
        <TileLayer url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png" />
        <GeoJSONLayer data={state.geojson as any} />
        <FitBounds bounds={state.bounds} />
        <InvalidateOnResize />
      </MapContainer>
    )
  );
};
