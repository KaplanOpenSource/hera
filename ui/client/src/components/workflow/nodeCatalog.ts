import { NodeCatalogEntry, WorkflowNode } from '../../shared/types';
import { FieldDef } from '../details/fieldDef';

// The entry type lives in shared, so the api package can use it too.
export type { NodeCatalogEntry };

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
export const isParameterKey = (name: string): boolean => {
  return /^[A-Za-z_][A-Za-z0-9_]*$/.test(name);
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

// The output names a node type produces, per the catalog. Empty when the type is
// unknown or has no detected outputs. Names only — the catalog carries no values.
export const nodeOutputNames = (
  node: WorkflowNode,
  catalog: NodeCatalogEntry[],
): string[] => {
  const entry = catalog.find(e => e.type === (node.type ?? ''));
  return (entry?.outputs ?? []).map(output => output.name);
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
      children[param.name] = { required: param.is_required, source: param.source };
    }
  }
  return { children };
};
