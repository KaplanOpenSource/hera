import L from 'leaflet';
import type { GeoTIFFImage } from 'geotiff';
import { tileRead, toGreyRgba } from './geoTiffTiles';

// Draws one band of a TIFF on a CRS.Simple map, reading only the blocks each tile needs.
export class GeoTiffGridLayer extends L.GridLayer {
  private image: GeoTIFFImage;
  private min: number;
  private max: number;

  constructor(image: GeoTIFFImage, min: number, max: number, options?: L.GridLayerOptions) {
    super(options);
    this.image = image;
    this.min = min;
    this.max = max;
  }

  createTile(coords: L.Coords, done: L.DoneCallback): HTMLElement {
    const tileSize = this.getTileSize();
    const canvas = document.createElement('canvas');
    canvas.width = tileSize.x;
    canvas.height = tileSize.y;
    const read = tileRead(coords.x, coords.y, coords.z, tileSize.x, this.image.getWidth(), this.image.getHeight());
    if (!read) {
      setTimeout(() => done(undefined, canvas), 0);
      return canvas;
    }
    this.drawTile(canvas, read.window, read.width, read.height).then(
      () => done(undefined, canvas),
      (error: Error) => {
        console.error('GeoTIFF tile read failed', error);
        done(error, canvas);
      },
    );
    return canvas;
  }

  private async drawTile(canvas: HTMLCanvasElement, window: number[], width: number, height: number) {
    const rasters = await this.image.readRasters({ window, width, height, samples: [0] });
    const band = rasters[0] as ArrayLike<number>;
    const rgba = toGreyRgba(band, this.min, this.max, this.image.getGDALNoData());
    canvas.getContext('2d')!.putImageData(new ImageData(rgba as Uint8ClampedArray<ArrayBuffer>, width, height), 0, 0);
  }
}
