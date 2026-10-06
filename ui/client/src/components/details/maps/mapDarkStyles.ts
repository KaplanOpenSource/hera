// OSM has no dark tiles, so dark mode inverts these instead (see darkMapControlsSx).
export const TILE_URL = 'https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png';
export const TILE_ATTRIBUTION = '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors';

// Darkens Leaflet's built-in zoom buttons and attribution box so they don't show as
// bright white boxes on a dark map. react-leaflet renders these into our own DOM
// (not a cross-origin iframe), so we can style them from here.
export const darkMapControlsSx = {
  // Leaflet's container defaults to a light gray, which flashes through while tiles
  // load and shows past the tile edges when zoomed out. Match the dark tiles.
  '& .leaflet-container': {
    backgroundColor: '#0e0e0e',
  },
  // Inverted, then the hue spun back, so water and parks keep roughly their colour.
  '& .leaflet-tile-pane': {
    filter: 'invert(1) hue-rotate(180deg) brightness(0.95) contrast(0.9)',
  },
  '& .leaflet-control-zoom a': {
    backgroundColor: '#2b2b2b',
    color: '#ddd',
    borderColor: '#000',
  },
  '& .leaflet-control-zoom a:hover': {
    backgroundColor: '#3a3a3a',
  },
  '& .leaflet-control-attribution': {
    backgroundColor: 'rgba(40, 40, 40, 0.85)',
    color: '#aaa',
  },
  '& .leaflet-control-attribution a': {
    color: '#8ab4f8',
  },
};
