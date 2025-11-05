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

export interface ProjectDocument {
  _cls: string;
  _id: { '$oid': string };
  projectName: string;
  desc: {
    toolkit?: string;
    datasourceName?: string;
    version?: number[];
    filesDirectory?: string;
    analysis_CacheCounter?: number;
    // docid?: string;
  },
  type: string;
  resource: string;
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
}