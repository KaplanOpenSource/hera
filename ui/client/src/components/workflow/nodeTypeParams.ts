import { isParameterKey, NodeCatalogEntry } from './nodeCatalog';

// A param counts as filled unless it's blank (null / undefined / ''). 0 and false
// are real values and are kept.
const hasValue = (value: any): boolean => {
  return value !== undefined && value !== null && value !== '';
};

// The input_parameters after a node's type changes to `entry`: keep every
// existing param the user gave a value, drop the empty ones, then seed the new
// type's params (empty) that aren't already kept. A param shared by the old and
// new type (e.g. ProjectName) keeps its value, since valued params are kept
// first. Existing params are returned unchanged for an unknown type.
export const paramsOnTypeChange = (
  existing: { [key: string]: any },
  entry: NodeCatalogEntry | undefined,
): { [key: string]: any } => {
  if (!entry) {
    return existing;
  }
  const next: { [key: string]: any } = {};
  for (const [key, value] of Object.entries(existing)) {
    if (hasValue(value)) {
      next[key] = value;
    }
  }
  for (const param of entry.parameters) {
    if (isParameterKey(param.name) && !(param.name in next)) {
      next[param.name] = '';
    }
  }
  return next;
};
