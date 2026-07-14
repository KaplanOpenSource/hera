// Generic, domain-agnostic definition of a field in the details tree. Recursive:
// `children` holds the defs of an object field's sub-fields, keyed by key — so
// one `def` describes a whole subtree. Extend with more details (type, default,
// help, …); callers pass one `def` and the field reads what it needs.
import { NodeParameterSource } from '../../shared/types';

export interface FieldDef {
  required?: boolean;
  // Where this field was discovered (e.g. a Hermes param's source). Shown as a
  // small badge on the row when present.
  source?: NodeParameterSource;
  children?: { [key: string]: FieldDef };
}
