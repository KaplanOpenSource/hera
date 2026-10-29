import { AgentConfig } from "./AgentConfig";

// Where the node catalog discovered a Hermes node's parameter or output.
export enum NodeParameterSource {
  JsonForm = 'jsonForm',
  Python = 'python',
  Template = 'template',
}

// One parameter a node type accepts. Values (defaults) are not provided by the
// catalog — only the name, whether it is required, and where it was discovered.
export interface NodeCatalogParameter {
  name: string;
  is_required: boolean;
  source: NodeParameterSource;
}

// One output a node type produces. Names only (no values), plus where it was
// discovered. Best-effort — only outputs detected from the node's run() return.
export interface NodeCatalogOutput {
  name: string;
  source: NodeParameterSource;
}

// One Hermes node type, the parameters it accepts and the outputs it produces.
export interface NodeCatalogEntry {
  type: string;
  parameters: NodeCatalogParameter[];
  // Optional: older callers/fixtures may omit it; the server always sends it.
  outputs?: NodeCatalogOutput[];
}

// The three collections a document can live in. The server names them alike.
export enum Collection {
  Measurements = 'Measurements',
  Simulations = 'Simulations',
  Cache = 'Cache',
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
  defaultTileServer?: string;
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
  // Flattened per-node input_parameters (keyed by node name). A query index the
  // database filters workflows on; kept in sync with the workflow by the UI.
  parameters?: { [node: string]: any };
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