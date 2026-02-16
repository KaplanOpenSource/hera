export type JsonValue = string | number | boolean | null | undefined | JsonValue[] | { [key: string]: JsonValue };

export interface CompareResult {
  path: string;
  values: (JsonValue | undefined)[];
}

function getAllPaths(obj: JsonValue, prefix = ""): string[] {
  if (obj === null || obj === undefined || typeof obj !== "object") {
    return [prefix || "/"];
  }
  const paths: string[] = [];
  if (Array.isArray(obj)) {
    paths.push(prefix || "/");
    obj.forEach((item, i) => paths.push(...getAllPaths(item, `${prefix}/${i}`)));
  } else {
    for (const key of Object.keys(obj)) {
      const path = `${prefix}/${key}`;
      if (typeof obj[key] === "object" && obj[key] !== null) {
        paths.push(...getAllPaths(obj[key], path));
      } else {
        paths.push(path);
      }
    }
  }
  return paths;
}

function getValueAtPath(obj: JsonValue, path: string): JsonValue | undefined {
  if (path === "/") return obj;
  const parts = path.split("/").filter(Boolean);
  let current: JsonValue = obj;
  for (const part of parts) {
    if (current === null || current === undefined || typeof current !== "object") return undefined;
    current = Array.isArray(current) ? current[parseInt(part, 10)] : (current as Record<string, JsonValue>)[part];
  }
  return current;
}

export function compareJsons(jsons: JsonValue[], omitEqual = false): CompareResult[] {
  const allPaths = new Set<string>();
  for (const json of jsons) {
    for (const p of getAllPaths(json)) allPaths.add(p);
  }

  const results: CompareResult[] = [];
  for (const path of allPaths) {
    const values = jsons.map((json) => getValueAtPath(json, path));

    if (omitEqual && values.every((v) => v === values[0])) {
      continue;
    }

    results.push({ path, values });
  }

  return results;
}

export function sortByDistinctValues(results: CompareResult[]): CompareResult[] {
  return [...results].sort((a, b) => {
    const distinctA = new Set(a.values).size;
    const distinctB = new Set(b.values).size;
    return distinctA - distinctB;
  });
}

export function filterAndSortByGroups(results: CompareResult[]): CompareResult[] {
  const withGroups = results.map((r) => {
    const nonUndef = r.values.filter((v) => v !== undefined);
    const groups = new Set(nonUndef).size;
    return { result: r, groups };
  });

  const filtered = withGroups.filter(({ groups }) => groups >= 2);
  const sorted = filtered.sort((a, b) => a.groups - b.groups);

  return sorted.map(({ result }) => result);
}
