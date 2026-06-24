import { WorkflowNode } from '../../shared/types';

// One parameter a node type accepts. Values (defaults) are not provided by the
// catalog — only the name, whether it is required, and where it was discovered.
export interface NodeCatalogParameter {
  name: string;
  is_required: boolean;
  source: string;
}

// One Hermes node type and the parameters it accepts.
export interface NodeCatalogEntry {
  type: string;
  parameters: NodeCatalogParameter[];
}

// The group a type belongs to in the picker: its dotted prefix without the last
// segment (e.g. `openFOAM.constant.g` → `openFOAM.constant`), or the whole type
// when it has no prefix (e.g. `BC`).
export const nodeTypeGroup = (type: string): string => {
  const dot = type.lastIndexOf('.');
  return dot === -1 ? type : type.slice(0, dot);
};

// A plain identifier we can use as an input_parameters key. The catalog also
// surfaces nested/template paths (e.g. `fields.items`) and loop variables that
// aren't real top-level parameters, so we keep only simple names.
const isParameterKey = (name: string): boolean => {
  return /^[A-Za-z_][A-Za-z0-9_]*$/.test(name);
};

// The input_parameters object to use after a node's type is set to `entry`.
// Each of the type's parameters is seeded so the user sees the shape to fill in;
// existing values are preserved (so re-picking a type never discards edited data)
// and the catalog gives names only, so new keys start empty.
export const prefilledParameters = (
  entry: NodeCatalogEntry | undefined,
  existing: { [key: string]: any },
): { [key: string]: any } => {
  if (!entry) {
    return existing;
  }
  const next: { [key: string]: any } = { ...existing };
  for (const param of entry.parameters) {
    if (isParameterKey(param.name) && !(param.name in next)) {
      next[param.name] = '';
    }
  }
  return next;
};

// Problems with a node, for a soft warning in the editor (never blocks editing).
// Returns an empty list when the node is fine, or while the catalog is still
// loading (we can't judge a type we don't have yet).
export const validateNode = (node: WorkflowNode, catalog: NodeCatalogEntry[]): string[] => {
  if (catalog.length === 0) {
    return [];
  }
  const type = node.type ?? '';
  if (!type) {
    return ['no type selected'];
  }
  const entry = catalog.find(e => e.type === type);
  if (!entry) {
    return [`unknown type: ${type}`];
  }
  const params = node.Execution?.input_parameters ?? {};
  const problems: string[] = [];
  for (const param of entry.parameters) {
    if (param.is_required && isParameterKey(param.name)) {
      const value = params[param.name];
      if (value === undefined || value === null || value === '') {
        problems.push(`missing required: ${param.name}`);
      }
    }
  }
  return problems;
};
