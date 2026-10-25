// Where a value sits inside a node's input_parameters: the keys from the top
// down to it, joined with a dot (e.g. `Parameters.project_name`, `Command.0`).
// A top-level parameter is a one segment path, which is its plain name.

// A list index is written as a path segment like any key.
const SEPARATOR = '.';

export const paramPathOf = (segments: string[]): string => {
  return segments.join(SEPARATOR);
};

export const segmentsOf = (path: string): string[] => {
  return path.split(SEPARATOR);
};

// The value at `path`, or undefined when any step of the way is missing.
export const valueAtPath = (params: { [key: string]: any }, path: string): any => {
  let current: any = params;
  for (const segment of segmentsOf(path)) {
    if (current === null || typeof current !== 'object') {
      return undefined;
    }
    current = current[segment];
  }
  return current;
};

// A shallow copy of `container`, still a list if it was one.
const copyOf = (container: any): any => {
  if (Array.isArray(container)) {
    return [...container];
  }
  if (container !== null && typeof container === 'object') {
    return { ...container };
  }
  return {};
};

// `params` with the value at `path` replaced. Every level on the way is copied,
// so the original is untouched. A missing level becomes a new object.
export const withValueAtPath = (
  params: { [key: string]: any },
  path: string,
  value: any,
): { [key: string]: any } => {
  const segments = segmentsOf(path);
  const root = copyOf(params);
  let container = root;
  for (const segment of segments.slice(0, -1)) {
    const next = copyOf(container[segment]);
    container[segment] = next;
    container = next;
  }
  container[segments[segments.length - 1]] = value;
  return root;
};

// `params` with the value at `path` gone: a key is deleted, a list element is
// spliced out. A path that is not there changes nothing.
export const withoutPath = (
  params: { [key: string]: any },
  path: string,
): { [key: string]: any } => {
  const segments = segmentsOf(path);
  const parentPath = paramPathOf(segments.slice(0, -1));
  const last = segments[segments.length - 1];
  if (parentPath === '') {
    const copy = { ...params };
    delete copy[last];
    return copy;
  }
  const parent = valueAtPath(params, parentPath);
  if (parent === null || typeof parent !== 'object') {
    return params;
  }
  if (Array.isArray(parent)) {
    return withValueAtPath(params, parentPath, parent.filter((_, i) => i !== Number(last)));
  }
  const copy = copyOf(parent);
  delete copy[last];
  return withValueAtPath(params, parentPath, copy);
};

// One string leaf and where it sits.
export interface ParamLeaf {
  path: string;
  text: string;
}

// Every string inside `params`, however deeply nested, with its path. Dicts and
// lists are both walked; other values (numbers, booleans, null) are skipped.
export const leafStringsOf = (params: { [key: string]: any }): ParamLeaf[] => {
  const found: ParamLeaf[] = [];
  const walk = (value: any, segments: string[]): void => {
    if (typeof value === 'string') {
      found.push({ path: paramPathOf(segments), text: value });
      return;
    }
    if (Array.isArray(value)) {
      value.forEach((item, index) => walk(item, [...segments, String(index)]));
      return;
    }
    if (value !== null && typeof value === 'object') {
      Object.entries(value).forEach(([key, item]) => walk(item, [...segments, key]));
    }
  };
  Object.entries(params).forEach(([key, value]) => walk(value, [key]));
  return found;
};
