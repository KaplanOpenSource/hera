// Generic, domain-agnostic definition of a field in the details tree. Recursive:
// `children` holds the defs of an object field's sub-fields, keyed by key — so
// one `def` describes a whole subtree. Extend with more details (type, default,
// help, …); callers pass one `def` and the field reads what it needs.
export interface FieldDef {
  required?: boolean;
  children?: { [key: string]: FieldDef };
}
