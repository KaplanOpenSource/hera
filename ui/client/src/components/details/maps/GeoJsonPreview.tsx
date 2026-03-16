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

const FitBounds = ({
  geojson,
}: {
  geojson: GeoJSON;
}) => {
  const map = useMap();
  useEffect(() => {
    const layer = L.geoJSON(geojson as any);
    const bounds = layer.getBounds();
    if (bounds.isValid()) {
      map.fitBounds(bounds, { padding: [20, 20] });
    }
  }, [geojson, map]);
  return null;
};

const loadGeoJson = async (path: string): Promise<GeoJSON> => {
  const { data, problem } = await fetchPython({
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
  if (problem) throw new Error(problem);
  return data.geojson_data;
};

export const GeoJsonPreview = ({
  path,
}: {
  path: string;
}) => {
  const [geojson, setGeojson] = useState<GeoJSON | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    loadGeoJson(path)
      .then(setGeojson)
      .catch(e => setError(e.message));
  }, [path]);

  if (error) {
    return <Typography color="error">Failed to load GeoJSON: {error}</Typography>;
  }
  if (!geojson) {
    return <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' }}><CircularProgress /></Box>;
  }

  return (
    <MapContainer
      center={[32, 35]}
      zoom={8}
      style={{ height: '100%', width: '100%' }}
    >
      <TileLayer url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png" />
      <GeoJSONLayer data={geojson as any} />
      <FitBounds geojson={geojson} />
      <InvalidateOnResize />
    </MapContainer>
  );
};
