import type { GeoJSON } from 'geojson';
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

// Every call to /init, /cors, /jupyter/ensure and /geojson. Each method
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

  // A GeoJSON file from the server's disk, converted to WGS84.
  async geoJson(path: string): Promise<GeoJSON | null> {
    return await this.getJson('/geojson', { path }, { label: 'geojson' });
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
