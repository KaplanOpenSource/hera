import { useEffect, useRef } from 'react';
import { useMap } from 'react-leaflet';

export const InvalidateOnResize = () => {
  const map = useMap();
  const boxRef = useRef(map.getContainer().parentElement);

  useEffect(() => {
    const el = boxRef.current;
    if (!el) return;
    const observer = new ResizeObserver(() => {
      map.invalidateSize();
    });
    observer.observe(el);
    return () => observer.disconnect();
  }, [map]);

  return null;
};
