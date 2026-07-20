import { Box, useTheme } from '@mui/material';
import { MapContainer, TileLayer } from 'react-leaflet';
import 'leaflet/dist/leaflet.css';
import { InvalidateOnResize } from './InvalidateOnResize';
import { darkMapControlsSx } from './mapDarkStyles';

export const isTileUrl = (resource: unknown): resource is string => {
  return typeof resource === 'string' && resource.includes('{x}') && resource.includes('{y}');
};

export const TileMapView = ({
  url,
}: {
  url: string,
}) => {
  // The tiles are the user's own data, so keep them; only darken the zoom and
  // attribution controls.
  const dark = useTheme().palette.mode === 'dark';
  return (
    <Box sx={{ height: '100%', width: '100%', ...(dark && darkMapControlsSx) }}>
      <MapContainer
        center={[32, 35]}
        zoom={8}
        style={{ height: '100%', width: '100%' }}
      >
        <TileLayer url={url} />
        <InvalidateOnResize />
      </MapContainer>
    </Box>
  );
};
