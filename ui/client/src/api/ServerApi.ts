import type { Feature, FeatureCollection } from 'geojson';
import { deserialize } from 'flatgeobuf/geojson';
import type { HeaderMeta } from 'flatgeobuf';
import { ProjectEntire, ProjectName } from '@shared/types';
import { ApiCall } from './ApiCall';

// What GET /init answers. Any piece may be missing: the server collects each
// one on its own and names the failures in `problems` (see init_data.py).
export type InitData = {
  ready: boolean,
  projects?: ProjectName[],
  datatypes?: { [key: string]: string },
  toolkitDocs?: any[],
  project?: ProjectEntire,
  server?: { heraVersion?: string, username?: string, inDocker?: boolean },
  problems?: { [piece: string]: string },
};

export type CorsInfo = {
  origins?: string[],
};

export type JupyterEnsureAnswer = {
  port: number,
};

// A map bbox in WGS84 degrees: west, south, east, north.
export type GeoBbox = [number, number, number, number];

// What GET /geoshapes answers, decoded. dotsLimit is set when the server sent
// one dot per shape: the most full shapes it sends for one bbox.
export type GeoShapesAnswer = {
  shapes: FeatureCollection,
  dotsLimit: number | null,
};

// The FlatGeobuf title the server writes on dots, see geo_shapes.dots_title.
const DOTS_TITLE = /^dots max=(\d+)$/;

// Every call to /init, /cors, /jupyter/ensure and /geoshapes. Each method
// answers data, or null when the call failed.
export class ServerApi extends ApiCall {
  // Everything the UI needs to start. Not ready yet carries no data.
  async init(wantedProjectName?: string): Promise<InitData | null> {
    const params: { [key: string]: string } = {};
    if (wantedProjectName) {
      params.project = wantedProjectName;
    }
    return await this.getJson('/init', params, { label: 'init' });
  }

  // A shapes file (.shp, .geojson, .fgb) from the server's disk, in WGS84.
  // With a bbox, only the full shapes in it, and no toasts: it runs on every map move.
  async geoShapes(path: string, bbox?: GeoBbox): Promise<GeoShapesAnswer | null> {
    const params: { [key: string]: string } = { path };
    let silent = false;
    if (bbox) {
      params.bbox = bbox.join(',');
      silent = true;
    }
    const bytes = await this.getBytes('/geoshapes', params, { label: 'geoshapes', silent });
    if (!bytes) {
      return null;
    }
    const features: Feature[] = [];
    // No bytes means no shapes: GDAL cannot write an empty FlatGeobuf.
    if (bytes.byteLength === 0) {
      return { shapes: { type: 'FeatureCollection', features }, dotsLimit: null };
    }
    let dotsLimit: number | null = null;
    const readTitle = (header: HeaderMeta) => {
      const match = header.title?.match(DOTS_TITLE);
      if (match) {
        dotsLimit = Number(match[1]);
      }
    };
    for await (const feature of deserialize(new Uint8Array(bytes), { headerMetaFn: readTitle })) {
      features.push(feature as Feature);
    }
    return { shapes: { type: 'FeatureCollection', features }, dotsLimit };
  }

  // The origins the server's CORS setup allows, for the header indicator.
  async cors(): Promise<CorsInfo | null> {
    return await this.getJson('/cors', {}, { label: 'cors', silent: true });
  }

  // Starts (or reuses) a Jupyter server for the given root directory.
  async ensureJupyter(rootDir: string, dark: boolean): Promise<JupyterEnsureAnswer | null> {
    return await this.postJson('/jupyter/ensure', { root_dir: rootDir, dark }, { label: 'jupyter ensure' });
  }
}
