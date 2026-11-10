// Pixel math and colors for drawing a TIFF as map tiles on a plain-pixel (CRS.Simple) map.

export type TileRead = {
  // Image pixels to read: [left, top, right, bottom].
  window: number[],
  // Canvas pixels to draw the read part into.
  width: number,
  height: number,
};

// The image area one tile covers, cut to the image; null when the tile is off the image.
export const tileRead = (
  x: number,
  y: number,
  zoom: number,
  tileSize: number,
  imageWidth: number,
  imageHeight: number,
): TileRead | null => {
  const span = tileSize / Math.pow(2, zoom);
  const left = Math.max(0, Math.floor(x * span));
  const top = Math.max(0, Math.floor(y * span));
  const right = Math.min(imageWidth, Math.ceil((x + 1) * span));
  const bottom = Math.min(imageHeight, Math.ceil((y + 1) * span));
  if (right <= left || bottom <= top) {
    return null;
  }
  return {
    window: [left, top, right, bottom],
    width: Math.max(1, Math.round((right - left) / span * tileSize)),
    height: Math.max(1, Math.round((bottom - top) / span * tileSize)),
  };
};

// Lowest and highest value, skipping no-data.
export const valueRange = (values: ArrayLike<number>, noData: number | null): { min: number, max: number } => {
  let min = Infinity;
  let max = -Infinity;
  for (let i = 0; i < values.length; i++) {
    const value = values[i];
    if (value === noData || Number.isNaN(value)) {
      continue;
    }
    if (value < min) {
      min = value;
    }
    if (value > max) {
      max = value;
    }
  }
  return { min, max };
};

// Grey from min (black) to max (white); no-data is see-through.
export const toGreyRgba = (
  values: ArrayLike<number>,
  min: number,
  max: number,
  noData: number | null,
): Uint8ClampedArray => {
  const rgba = new Uint8ClampedArray(values.length * 4);
  let scale = 0;
  if (max > min) {
    scale = 255 / (max - min);
  }
  for (let i = 0; i < values.length; i++) {
    const value = values[i];
    if (value === noData || Number.isNaN(value)) {
      continue;
    }
    const grey = (value - min) * scale;
    rgba[i * 4] = grey;
    rgba[i * 4 + 1] = grey;
    rgba[i * 4 + 2] = grey;
    rgba[i * 4 + 3] = 255;
  }
  return rgba;
};
