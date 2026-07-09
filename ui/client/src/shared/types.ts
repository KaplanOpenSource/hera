import { AgentConfig } from "./AgentConfig";

// Where the node catalog discovered a Hermes node's parameter or output.
export enum NodeParameterSource {
  JsonForm = 'jsonForm',
  Python = 'python',
  Template = 'template',
}

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

export interface DocumentDesc {
  toolkit?: string;
  datasourceName?: string;
  version?: number[] | string;
  filesDirectory?: string;
  analysis_CacheCounter?: number;
  // docid?: string;
}

export interface WorkflowNode {
  type?: string;
  // Names of node(s) that must run before this one. Hermes accepts a single
  // name or a list.
  requires?: string | string[];
  Execution?: { input_parameters?: { [key: string]: any } };
}

export interface WorkflowBlock {
  root?: any;
  solver?: string;
  nodeList?: string[];
  nodes?: { [name: string]: WorkflowNode };
  [key: string]: any;
}

// desc.workflow holds the block, optionally wrapped in one extra { workflow }
// level (the shape Hermes writes).
export type WorkflowData = WorkflowBlock | { workflow?: WorkflowBlock };

// Desc of a Hermes workflow document — the base desc plus workflow-specific
// fields, kept off DocumentDesc so the base doesn't accumulate optionals.
export interface WorkflowDesc extends DocumentDesc {
  workflowName?: string;
  workflow?: WorkflowData;
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

export interface Problem {
  error: string;
  traceback: string;
}

export interface ExecResponse {
  data: any;
  problem: Problem | null;
}

export interface Repository {
  datasourceName: string;
  resource: string;
  dataFormat: string;
  toolkit: string;
  version: number[];
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