import { describe, it, expect } from 'vitest';
import { tileRead, toGreyRgba, valueRange } from '../src/components/details/maps/geoTiffTiles';

describe('tileRead', () => {
  it('at zoom 0 a tile is 256 image pixels', () => {
    expect(tileRead(1, 2, 0, 256, 1000, 1000)).toEqual({ window: [256, 512, 512, 768], width: 256, height: 256 });
  });

  it('zoomed out, one tile covers more image pixels', () => {
    expect(tileRead(0, 0, -2, 256, 5000, 5000)).toEqual({ window: [0, 0, 1024, 1024], width: 256, height: 256 });
  });

  it('cuts a tile at the image edge and shrinks the canvas part to match', () => {
    expect(tileRead(0, 0, -2, 256, 512, 2000)).toEqual({ window: [0, 0, 512, 1024], width: 128, height: 256 });
  });

  it('is null for a tile off the image', () => {
    expect(tileRead(5, 0, 0, 256, 1000, 1000)).toBeNull();
    expect(tileRead(-1, 0, 0, 256, 1000, 1000)).toBeNull();
  });
});

describe('valueRange', () => {
  it('skips no-data', () => {
    expect(valueRange([3, 9, 255, 1], 255)).toEqual({ min: 1, max: 9 });
  });
});

describe('toGreyRgba', () => {
  it('maps min to black, max to white, no-data to see-through', () => {
    const rgba = toGreyRgba([0, 16, 99], 0, 16, 99);

    expect([...rgba]).toEqual([0, 0, 0, 255, 255, 255, 255, 255, 0, 0, 0, 0]);
  });

  it('draws a single value image black, not blank', () => {
    expect([...toGreyRgba([5], 5, 5, null)]).toEqual([0, 0, 0, 255]);
  });
});
