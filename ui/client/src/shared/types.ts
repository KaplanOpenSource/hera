export interface ProjectName {
  id: string;
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

// The /exec endpoint returns a plain JSON value
export interface ExecRequest {
  code: string;
}


