import { Box } from '@mui/material';
import { MapContainer, TileLayer } from 'react-leaflet';
import 'leaflet/dist/leaflet.css';

export const isTileUrl = (resource: unknown): resource is string => {
  return typeof resource === 'string' && resource.includes('{x}') && resource.includes('{y}');
};

export const TileMapView = ({
  url,
}: {
  url: string,
}) => {
  return (
    <Box sx={{ flex: 1, minHeight: 0, mx: -2, mb: -2 }}>
      <MapContainer
        center={[32, 35]}
        zoom={8}
        style={{ height: '100%', width: '100%' }}
      >
        <TileLayer url={url} />
      </MapContainer>
    </Box>
  );
};
