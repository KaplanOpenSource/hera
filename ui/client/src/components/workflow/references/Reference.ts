import { ReferenceKind } from './ReferenceKind';

// Prefix on dataflow edge ids: df:<node>:<handleMark>:<key>-><target>.<param>.
export const DATAFLOW_EDGE_PREFIX = 'df:';

// One reference, as a value: whose it is, of what kind, and which key.
export class Reference {
  constructor(
    readonly node: string,
    readonly kind: ReferenceKind,
    readonly key: string,
  ) {}

  toString(): string {
    return this.kind.format(this.node, this.key);
  }

  handleId(): string {
    return this.kind.handleId(this.node, this.key);
  }

  // The id of the canvas line from this reference to the parameter reading it.
  edgeIdTo(target: string, param: string): string {
    return `${DATAFLOW_EDGE_PREFIX}${this.handleId()}->${target}.${param}`;
  }

  clearToken(): RegExp {
    return this.kind.clearToken(this.node, this.key);
  }
}
