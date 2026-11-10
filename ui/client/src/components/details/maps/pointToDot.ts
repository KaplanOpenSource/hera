import L from 'leaflet';

const DOT_RADIUS = 3;

// Small canvas circles, not Leaflet's pin icons, so many points stay usable.
export const pointToDot = (_feature: unknown, latlng: L.LatLng): L.Layer => {
  return L.circleMarker(latlng, { radius: DOT_RADIUS });
};
