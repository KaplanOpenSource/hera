import { AgentConfig } from "./AgentConfig";

export interface ProjectName {
  name: string;
}

export interface Project {
  id: string;
  name: string;
  documents?: {
    sim: number;
    measure: number;
    cache: number;
  };
  toolkitCount?: number;
}

export interface ProjectEntire {
  name: string;
  documents: ProjectDocument[];
}

export enum MetadataCls {
  Simulations = 'Metadata.Simulations',
  Measurements = 'Metadata.Measurements',
  Cache = 'Metadata.Cache',
}

export interface DocumentDesc {
  toolkit?: string;
  datasourceName?: string;
  version?: number[] | string;
  filesDirectory?: string;
  analysis_CacheCounter?: number;
  // docid?: string;
}

export interface ProjectDocument {
  _cls: string;
  _id: { '$oid': string };
  projectName: string;
  desc: DocumentDesc,
  type: string;
  resource: string | AgentConfig;
  dataFormat: string;
}

// The /exec endpoint returns a plain JSON value
export interface ExecRequest {
  code: string;
}

export interface Toolkit {
  toolkit: string;
  cls: string;
  description?: string;
  source?: string;
  type?: string;
  repositoryName?: string;
  version?: string;
}