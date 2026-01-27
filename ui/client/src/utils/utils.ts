/**
 * Copies an object without the specified fields
 * @param obj Object to be cloned
 * @param fields Fields to omit from output
 * @returns The cloned object without the fields
 */
export const copyWithout = (obj: any, fields: string[]) => {
  return Object.fromEntries(Object.entries(obj).filter(([k, _v]) => !fields.includes(k)));
};
