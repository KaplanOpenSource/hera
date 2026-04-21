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
    label: 'load GeoJSON',
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
  const [mapState, setMapState] = useState<MapState>({ geojson: null, bounds: null, hasError: false });

  useEffect(() => {
    (async () => {
      setMapState({ geojson: null, bounds: null, hasError: false });
      const geojson = await loadGeoJson(path);
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
      <MapContainer
        center={[32, 35]}
        zoom={8}
        style={{ height: '100%', width: '100%' }}
      >
        <TileLayer url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png" />
        <GeoJSONLayer data={mapState.geojson as any} />
        <FitBounds bounds={mapState.bounds} />
        <InvalidateOnResize />
      </MapContainer>
    )
  );
};
