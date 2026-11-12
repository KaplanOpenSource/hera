import { ReferenceKind } from './ReferenceKind';

// Prefix on dataflow edge ids: df:<node>:<handleMark>:<key>-><target>.<path>.
export const DATAFLOW_EDGE_PREFIX = 'df:';

// One reference, as a value: whose it is, of what kind, and which key.
export class Reference {
  constructor(
    readonly node: string,
    readonly kind: ReferenceKind,
    readonly key: string,
  ) {}

  // The output's own name, without any JSONPath into it (`items[0].name` -> `items`).
  rootKey(): string {
    return this.key.split(/[.[]/)[0];
  }

  toString(): string {
    return this.kind.format(this.node, this.key);
  }

  handleId(): string {
    return this.kind.handleId(this.node, this.key);
  }

  // The dot a line actually leaves from: one per output, shared by its sub-fields.
  rootHandleId(): string {
    return this.kind.handleId(this.node, this.rootKey());
  }

  // The id of the canvas line from this reference to the parameter reading it.
  edgeIdTo(target: string, paramPath: string): string {
    return `${DATAFLOW_EDGE_PREFIX}${this.handleId()}->${target}.${paramPath}`;
  }

  clearToken(): RegExp {
    return this.kind.clearToken(this.node, this.key);
  }
}
