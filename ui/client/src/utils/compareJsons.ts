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

export function getValueAtPath(obj: JsonValue, path: string): JsonValue | undefined {
  if (path === "/") return obj;
  const parts = path.split("/").filter(Boolean);
  let current: JsonValue = obj;
  for (const part of parts) {
    if (current === null || current === undefined || typeof current !== "object") return undefined;
    current = Array.isArray(current) ? current[parseInt(part, 10)] : (current as Record<string, JsonValue>)[part];
  }
  return current;
}

export function compareJsons(jsonsRaw: any[], omitEqual = false): CompareResult[] {
  const jsons: JsonValue[] = jsonsRaw;
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

export function filterAndSortByGroups(results: CompareResult[], minGroupSize: number): CompareResult[] {
  const withGroups: { result: CompareResult; distinctValues: number }[] = [];

  for (const r of results) {
    const valueCounts = new Map<JsonValue, number>();
    for (const v of r.values) {
      if (v === undefined) continue;
      valueCounts.set(v, (valueCounts.get(v) || 0) + 1);
    }
    const maxCount = Math.max(0, ...valueCounts.values());
    if (valueCounts.size >= 2 && maxCount >= minGroupSize) {
      withGroups.push({ result: r, distinctValues: valueCounts.size });
    }
  }

  const sorted = withGroups.sort((a, b) => a.distinctValues - b.distinctValues);
  return sorted.map(({ result }) => result);
}
