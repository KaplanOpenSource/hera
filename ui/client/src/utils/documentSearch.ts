// Collect all primitive values from a document into one lowercased string for substring
// search. Keys and JSON punctuation are intentionally excluded, so matches are on values only.
const collectValues = (value: unknown, out: string[]): void => {
  if (value === null || value === undefined) return;
  if (Array.isArray(value)) {
    for (const item of value) collectValues(item, out);
  } else if (typeof value === 'object') {
    for (const item of Object.values(value as { [key: string]: unknown })) collectValues(item, out);
  } else {
    out.push(String(value));
  }
};

export const documentSearchText = (data: unknown): string => {
  const out: string[] = [];
  collectValues(data, out);
  return out.join(' ').toLowerCase();
};

// The document fields a scoped `field:value` term can target. The text before the
// first colon must start with one of these to be treated as a field query.
const FIELD_ROOTS = ['type', 'resource', 'desc'];

// A term's left side is only treated as a field name when it looks like a dotted
// identifier path (e.g. `type`, `desc.groupName`), so that values containing a
// colon such as `12:30` stay plain-text searches.
const IDENTIFIER_PATH = /^[A-Za-z_][A-Za-z0-9_]*(\.[A-Za-z_][A-Za-z0-9_]*)*$/;

export enum SearchTermKind {
  Plain = 'plain',
  Field = 'field',
  Unknown = 'unknown',
}

export interface SearchTerm {
  kind: SearchTermKind;
  value: string;      // lowercased text to match
  path?: string[];    // resolved field path, for Field terms
  field?: string;     // the (unknown) root name, for Unknown terms
}

// Splits a raw query into whitespace-separated terms, classifying each as a plain
// text search, a scoped field search, or a search on an unknown field.
export const parseSearchQuery = (query: string): SearchTerm[] => {
  const parseTerm = (raw: string): SearchTerm => {
    const colon = raw.indexOf(':');
    if (colon > 0) {
      const left = raw.slice(0, colon);
      const value = raw.slice(colon + 1).toLowerCase();
      if (IDENTIFIER_PATH.test(left)) {
        const segments = left.split('.');
        const root = segments[0].toLowerCase();
        if (FIELD_ROOTS.includes(root)) {
          return { kind: SearchTermKind.Field, value, path: [root, ...segments.slice(1)] };
        }
        return { kind: SearchTermKind.Unknown, value, field: segments[0] };
      }
    }
    return { kind: SearchTermKind.Plain, value: raw.toLowerCase() };
  };
  return query.trim().split(/\s+/).filter(Boolean).map(parseTerm);
};

// The unknown field names referenced by a parsed query (deduplicated), for warning the user.
export const unknownSearchFields = (terms: SearchTerm[]): string[] => {
  const names = terms.filter(t => t.kind === SearchTermKind.Unknown).map(t => t.field!);
  return [...new Set(names)];
};

// Resolves a dotted field path against the document data.
const valueAtPath = (data: unknown, path: string[]): unknown => {
  let current: unknown = data;
  for (const segment of path) {
    if (current === null || current === undefined || typeof current !== 'object') return undefined;
    current = (current as { [key: string]: unknown })[segment];
  }
  return current;
};

// The lowercased value blob of a single subtree, for scoped field matching.
const subtreeText = (value: unknown): string => {
  const out: string[] = [];
  collectValues(value, out);
  return out.join(' ').toLowerCase();
};

// Whether a document matches every term (AND). `blob` is the precomputed full-document
// value text used for plain terms; field terms are resolved against `data` directly.
export const documentMatchesQuery = (data: unknown, blob: string, terms: SearchTerm[]): boolean => {
  const matchesTerm = (term: SearchTerm): boolean => {
    if (term.kind === SearchTermKind.Plain) return blob.includes(term.value);
    if (term.kind === SearchTermKind.Unknown) return false;
    const resolved = valueAtPath(data, term.path!);
    return resolved === undefined ? false : subtreeText(resolved).includes(term.value);
  };
  return terms.every(matchesTerm);
};
