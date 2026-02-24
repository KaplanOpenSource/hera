export type JsonValue = string | number | boolean | null | undefined | JsonValue[] | { [key: string]: JsonValue };

export interface CompareResult {
  path: string;
  values: (JsonValue | undefined)[];
}

/**
 * Recursively collects all leaf paths in a JSON value using slash-delimited format (e.g., "/foo/bar/0").
 * For primitives, returns ["/"]. For arrays, includes the array path itself plus its children.
 * For objects, only includes paths to leaf (non-object) values.
 */
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

/**
 * Retrieves the value at a slash-delimited path (e.g., "/foo/0/bar").
 * Handles both object keys and array indices. Returns undefined if
 * the path doesn't exist or traverses through a primitive.
 */
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

/**
 * Compares multiple JSON values by collecting all unique paths across them
 * and returning each path with its corresponding value from every input.
 * If omitEqual is true, paths where all values are identical are excluded.
 *
 * @example
 * compareJsons(
 *   [
 *     { name: "alice", age: 30, city: "NYC" },
 *     { name: "bob",   age: 30, city: "LA" },
 *     { name: "carol", age: 30, city: "NYC" },
 *   ],
 *   true
 * )
 * // => [
 * //   { path: "/name", values: ["alice", "bob", "carol"] },
 * //   { path: "/city", values: ["NYC", "LA", "NYC"] },
 * // ]
 * // "/age" is omitted because all three values are 30
 */
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

/**
 * Sorts comparison results by the number of distinct values ascending,
 * surfacing paths with the most agreement first.
 * Note: uses reference equality via Set, so objects/arrays won't deduplicate as expected.
 */
export function sortByDistinctValues(results: CompareResult[]): CompareResult[] {
  return [...results].sort((a, b) => {
    const distinctA = new Set(a.values).size;
    const distinctB = new Set(b.values).size;
    return distinctA - distinctB;
  });
}

/**
 * Filters results to paths that have at least 2 distinct defined values where
 * the most common value appears at least minGroupSize times and there are
 * no more than maxBranches distinct values, then sorts by number of distinct
 * values ascending. Useful for finding paths that cleanly partition the inputs
 * into a small number of groups.
 *
 * @example
 * filterAndSortByGroups(
 *   [
 *     { path: "/name", values: ["alice", "bob", "carol", "dave"] },  // 4 distinct, max group size 1
 *     { path: "/city", values: ["NYC", "NYC", "LA", "LA"] },         // 2 distinct, max group size 2
 *     { path: "/role", values: ["admin", "admin", "admin", "user"] }, // 2 distinct, max group size 3
 *   ],
 *   2,
 *   3
 * )
 * // => [
 * //   { path: "/city", values: ["NYC", "NYC", "LA", "LA"] },
 * //   { path: "/role", values: ["admin", "admin", "admin", "user"] },
 * // ]
 * // "/name" is excluded because its largest group has size 1 (< minGroupSize 2)
 * // and it has 4 distinct values (> maxBranches 3)
 */
export function filterAndSortByGroups(results: CompareResult[], minGroupSize: number, maxBranches?: number): CompareResult[] {
  const withGroups: { result: CompareResult; distinctValues: number }[] = [];

  for (const r of results) {
    const valueCounts = new Map<JsonValue, number>();
    for (const v of r.values) {
      if (v === undefined) continue;
      valueCounts.set(v, (valueCounts.get(v) || 0) + 1);
    }
    const maxCount = Math.max(0, ...valueCounts.values());
    if (valueCounts.size >= 2 && maxCount >= minGroupSize) {
      if (maxBranches !== undefined && valueCounts.size > maxBranches) continue;
      withGroups.push({ result: r, distinctValues: valueCounts.size });
    }
  }

  const sorted = withGroups.sort((a, b) => a.distinctValues - b.distinctValues);
  return sorted.map(({ result }) => result);
}
