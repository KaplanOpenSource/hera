// Collect all primitive values from a document into one lowercased string for substring
// search. Keys and JSON punctuation are intentionally excluded, so matches are on values only.
const collectValues = (value: unknown, out: string[]): void => {
  if (value === null || value === undefined) return;
  if (Array.isArray(value)) {
    for (const item of value) collectValues(item, out);
  } else if (typeof value === 'object') {
    for (const item of Object.values(value as Record<string, unknown>)) collectValues(item, out);
  } else {
    out.push(String(value));
  }
};

export const documentSearchText = (data: unknown): string => {
  const out: string[] = [];
  collectValues(data, out);
  return out.join(' ').toLowerCase();
};
