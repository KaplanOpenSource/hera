/**
 * Copies an object without the specified fields
 * @param obj Object to be cloned
 * @param fields Fields to omit from output
 * @returns The cloned object without the fields
 */
export const copyWithout = (obj: any, fields: string[]) => {
  return Object.fromEntries(Object.entries(obj).filter(([k, _v]) => !fields.includes(k)));
};

/**
 * Reorder the output of `Object.entries()`
 * @param entries output of `Object.entries()` as array of key value
 * @param firstFields fields to move to the front
 * @param lastFields fields to move to the back
 * @returns reordered entries
 */
export const reorderEntries = (entries: [string, unknown][], firstFields: string[] = [], lastFields: string[] = []) => {
  return [
    ...entries.filter(([k]) => firstFields.includes(k)),
    ...entries.filter(([k]) => !firstFields.includes(k) && !lastFields.includes(k)),
    ...entries.filter(([k]) => lastFields.includes(k)),
  ];
};

/**
 * Returns a new array reordered so that items matching orderedKeys (by keySelector)
 * appear first in the specified order, followed by all remaining items in their
 * original order.
 *
 * @example
 * reorderByKeys(
 *   [{ id: "c" }, { id: "a" }, { id: "b" }, { id: "d" }],
 *   (x) => x.id,
 *   ["b", "a"]
 * )
 * // => [{ id: "b" }, { id: "a" }, { id: "c" }, { id: "d" }]
 */
export function reorderByKeys<T, K>(
  array: T[],
  keySelector: (item: T) => K,
  orderedKeys: K[]
): T[] {
  const keyIndex = new Map(orderedKeys.map((k, i) => [k, i]));
  const front: T[][] = orderedKeys.map(() => []);
  const rest: T[] = [];

  for (const item of array) {
    const idx = keyIndex.get(keySelector(item));
    if (idx !== undefined) {
      front[idx].push(item);
    } else {
      rest.push(item);
    }
  }

  return front.flat().concat(rest);
}