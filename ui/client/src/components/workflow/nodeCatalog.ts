import { WorkflowNode } from '../../shared/types';
import { FieldDef } from '../details/fieldDef';

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

// The type-level problem with a node, for a soft warning (never blocks editing).
// undefined when the type is fine, or while the catalog is still loading.
export const nodeTypeIssue = (node: WorkflowNode, catalog: NodeCatalogEntry[]): string | undefined => {
  if (catalog.length === 0) {
    return undefined;
  }
  const type = node.type ?? '';
  if (!type) {
    return 'no type selected';
  }
  if (!catalog.some(entry => entry.type === type)) {
    return `unknown type: ${type}`;
  }
  return undefined;
};

// The field def for a node's input_parameters: a `children` map (keyed by
// parameter name) marking which params are required. The single place that maps
// the Hermes catalog to the generic FieldDef — extend the mapping here as the
// catalog gains detail (type, default, …). Skips the noisy nested/template
// names, same as prefill.
export const paramsFieldDef = (
  node: WorkflowNode,
  catalog: NodeCatalogEntry[],
): FieldDef => {
  const entry = catalog.find(e => e.type === (node.type ?? ''));
  const children: { [key: string]: FieldDef } = {};
  for (const param of entry?.parameters ?? []) {
    if (isParameterKey(param.name)) {
      children[param.name] = { required: param.is_required };
    }
  }
  return { children };
};
