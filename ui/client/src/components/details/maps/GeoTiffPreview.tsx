import { useEffect, useState } from 'react';
import { Box, CircularProgress, Typography, useTheme } from '@mui/material';
import { MapContainer, useMap } from 'react-leaflet';
import L from 'leaflet';
import { fromUrl, type GeoTIFFImage } from 'geotiff';
import 'leaflet/dist/leaflet.css';
import { BASEURL } from '../../../shared/baseurl';
import { InvalidateOnResize } from './InvalidateOnResize';
import { GeoTiffGridLayer } from './GeoTiffGridLayer';
import { valueRange } from './geoTiffTiles';
import { darkMapControlsSx } from './mapDarkStyles';

const TIFF_EXTENSIONS = /\.(tif|tiff)$/i;

// Long side of the coarse read used to find the value range.
const RANGE_SAMPLE_SIZE = 512;

// Zoom limits for the map and its layer; Leaflet layers draw nothing below zoom 0 unless told.
const MIN_ZOOM = -10;
const MAX_ZOOM = 4;

export const isGeoTiffResource = (resource: unknown): resource is string => {
  return typeof resource === 'string' && TIFF_EXTENSIONS.test(resource);
};

type TiffState = {
  image: GeoTIFFImage | null,
  min: number,
  max: number,
  hasError: boolean,
};

const EMPTY_STATE: TiffState = { image: null, min: 0, max: 0, hasError: false };

// Pixel bounds of the image on a CRS.Simple map: row r sits at lat -r.
const imageBounds = (image: GeoTIFFImage): L.LatLngBounds => {
  return L.latLngBounds([[-image.getHeight(), 0], [0, image.getWidth()]]);
};

const TiffLayer = ({
  image,
  min,
  max,
}: {
  image: GeoTIFFImage,
  min: number,
  max: number,
}) => {
  const map = useMap();
  useEffect(() => {
    const bounds = imageBounds(image);
    const layer = new GeoTiffGridLayer(image, min, max, { bounds, noWrap: true, minZoom: MIN_ZOOM, maxZoom: MAX_ZOOM });
    layer.addTo(map);
    map.fitBounds(bounds);
    return () => {
      layer.remove();
    };
  }, [image, min, max, map]);
  return null;
};

export const GeoTiffPreview = ({
  path,
}: {
  path: string,
}) => {
  const [tiffState, setTiffState] = useState<TiffState>(EMPTY_STATE);
  const dark = useTheme().palette.mode === 'dark';

  useEffect(() => {
    let stale = false;
    setTiffState(EMPTY_STATE);
    const url = `${BASEURL}/file${path.startsWith('/') ? '' : '/'}${path}`;
    (async () => {
      const tiff = await fromUrl(url);
      const image = await tiff.getImage();
      const scale = RANGE_SAMPLE_SIZE / Math.max(image.getWidth(), image.getHeight());
      const sample = await image.readRasters({
        samples: [0],
        width: Math.max(1, Math.round(image.getWidth() * scale)),
        height: Math.max(1, Math.round(image.getHeight() * scale)),
      });
      const { min, max } = valueRange(sample[0] as ArrayLike<number>, image.getGDALNoData());
      if (!stale) {
        setTiffState({ image, min, max, hasError: false });
      }
    })().catch((error: Error) => {
      // A bad or missing file shows a message instead of a spinner that never ends.
      console.error('GeoTIFF load failed', error);
      if (!stale) {
        setTiffState({ ...EMPTY_STATE, hasError: true });
      }
    });
    return () => {
      stale = true;
    };
  }, [path]);

  const centered = { display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' };

  return (
    tiffState.hasError ? (
      <Box sx={centered}><Typography color="text.secondary">TIFF unavailable</Typography></Box>
    ) : !tiffState.image ? (
      <Box sx={centered}><CircularProgress /></Box>
    ) : (
      // The TIFF is data, so its tiles are never inverted in dark mode.
      <Box sx={{ height: '100%', width: '100%', ...(dark && darkMapControlsSx), '& .leaflet-tile-pane': { filter: 'none' } }}>
        <MapContainer
          crs={L.CRS.Simple}
          bounds={imageBounds(tiffState.image)}
          minZoom={MIN_ZOOM}
          maxZoom={MAX_ZOOM}
          style={{ height: '100%', width: '100%' }}
        >
          <TiffLayer image={tiffState.image} min={tiffState.min} max={tiffState.max} />
          <InvalidateOnResize />
        </MapContainer>
      </Box>
    )
  );
};
