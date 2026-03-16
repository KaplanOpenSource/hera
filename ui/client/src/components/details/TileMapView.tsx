import { Box } from '@mui/material';
import { useEffect, useRef } from 'react';
import { MapContainer, TileLayer } from 'react-leaflet';
import type { Map } from 'leaflet';
import 'leaflet/dist/leaflet.css';

export const isTileUrl = (resource: unknown): resource is string => {
  return typeof resource === 'string' && resource.includes('{x}') && resource.includes('{y}');
};

export const TileMapView = ({
  url,
}: {
  url: string,
}) => {
  const mapRef = useRef<Map | null>(null);
  const boxRef = useRef<HTMLDivElement | null>(null);

  useEffect(() => {
    if (!boxRef.current || !mapRef.current) return;
    const observer = new ResizeObserver(() => {
      mapRef.current?.invalidateSize();
    });
    observer.observe(boxRef.current);
    return () => observer.disconnect();
  }, []);

  return (
    <Box ref={boxRef} sx={{ height: '100%' }}>
      <MapContainer
        ref={mapRef}
        center={[32, 35]}
        zoom={8}
        style={{ height: '100%', width: '100%' }}
      >
        <TileLayer url={url} />
      </MapContainer>
    </Box>
  );
};
