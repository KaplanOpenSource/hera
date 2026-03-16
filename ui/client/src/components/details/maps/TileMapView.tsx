import { VisibilityOff } from '@mui/icons-material';
import { Box } from '@mui/material';
import { MapContainer, TileLayer } from 'react-leaflet';
import 'leaflet/dist/leaflet.css';
import { ButtonTooltip } from '../../../elements/ButtonTooltip';
import { InvalidateOnResize } from './InvalidateOnResize';

export const isTileUrl = (resource: unknown): resource is string => {
  return typeof resource === 'string' && resource.includes('{x}') && resource.includes('{y}');
};

export const TileMapView = ({
  url,
  onClose,
}: {
  url: string,
  onClose: () => void,
}) => {
  return (
    <Box sx={{ height: '100%', position: 'relative' }}>
      <Box sx={{ position: 'absolute', top: 4, right: 4, zIndex: 1000 }}>
        <ButtonTooltip
          title="Hide preview"
          onClick={onClose}
          sx={{
            backgroundColor: 'white',
            '&:hover': { backgroundColor: '#eee' },
          }}
        >
          <VisibilityOff sx={{ fontSize: 14 }} />
        </ButtonTooltip>
      </Box>
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
