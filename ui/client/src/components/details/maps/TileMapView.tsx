import { Box } from '@mui/material';
import { MapContainer, TileLayer } from 'react-leaflet';
import 'leaflet/dist/leaflet.css';
import { InvalidateOnResize } from './InvalidateOnResize';

export const isTileUrl = (resource: unknown): resource is string => {
  return typeof resource === 'string' && resource.includes('{x}') && resource.includes('{y}');
};

export const TileMapView = ({
  url,
}: {
  url: string,
}) => {
  return (
    <MapContainer
      center={[32, 35]}
      zoom={8}
      style={{ height: '100%', width: '100%' }}
    >
      <TileLayer url={url} />
      <InvalidateOnResize />
    </MapContainer>
  );
};
