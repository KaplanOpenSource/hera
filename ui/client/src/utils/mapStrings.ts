// A copy of a value with every string in it rewritten. Keys are left alone.
export const mapStrings = (value: any, rewrite: (text: string) => string): any => {
  if (typeof value === 'string') {
    return rewrite(value);
  }
  if (Array.isArray(value)) {
    return value.map(item => mapStrings(item, rewrite));
  }
  if (typeof value === 'object' && value !== null) {
    const mapped: { [key: string]: any } = {};
    for (const key of Object.keys(value)) {
      mapped[key] = mapStrings(value[key], rewrite);
    }
    return mapped;
  }
  return value;
};
