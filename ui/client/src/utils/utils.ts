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

