import { WorkflowNode } from '../../../shared/types';
import { NodeCatalogEntry } from '../nodeCatalog';
import { Reference } from './Reference';
import { ReferenceKind } from './ReferenceKind';

// Any reference in a parameter value. The section is greedy so a kind whose
// section holds a dot still parses; the last dotted part is always the key.
const PARSE = /\{\s*(\w+)\.([\w.]+)\.(\w+)\s*\}/g;
// A source dot's handle id: <node>:<handleMark>:<key>.
const HANDLE = /^(\w+):(\w+):(.+)$/;
// A dataflow edge id: df:<node>:<handleMark>:<key>-><target>.<paramPath>.
// The path may itself hold dots, so it is everything after the target's dot.
const EDGE_ID = /^df:(\w+):(\w+):(\w+)->(\w+)\.(.+)$/;

// The questions that are about all the kinds at once rather than one of them.
export class ReferenceKindRegistry {
  private readonly kinds: ReferenceKind[];
  private readonly byMark: { [mark: string]: ReferenceKind } = {};

  constructor(kinds: ReferenceKind[]) {
    this.kinds = kinds;
    for (const kind of kinds) {
      this.byMark[kind.handleMark] = kind;
    }
  }

  all(): ReferenceKind[] {
    return this.kinds;
  }

  // The kind a written section belongs to, or null when no kind reads it.
  bySection(section: string): ReferenceKind | null {
    for (const kind of this.kinds) {
      if (new RegExp(`^(?:${kind.sectionMatch()})$`).test(section)) {
        return kind;
      }
    }
    return null;
  }

  // Every reference written in one parameter value. Unknown sections are skipped.
  parseAll(value: string): Reference[] {
    const found: Reference[] = [];
    PARSE.lastIndex = 0;
    let match = PARSE.exec(value);
    while (match !== null) {
      const kind = this.bySection(match[2]);
      if (kind !== null) {
        found.push(new Reference(match[1], kind, match[3]));
      }
      match = PARSE.exec(value);
    }
    return found;
  }

  // One value's text with every known reference to `oldName` pointing at `newName`.
  renamedNode(value: string, oldName: string, newName: string): string {
    return value.replace(PARSE, (match, node, section, key) => {
      if (node !== oldName || this.bySection(section) === null) {
        return match;
      }
      return `{${newName}.${section}.${key}}`;
    });
  }

  // The reference a source dot stands for, or null when the id is not one.
  ofHandle(handleId: string): Reference | null {
    const match = HANDLE.exec(handleId);
    if (match === null) {
      return null;
    }
    const kind = this.byMark[match[2]];
    if (kind === undefined) {
      return null;
    }
    return new Reference(match[1], kind, match[3]);
  }

  // The reference and the parameter path a dataflow edge id links, or null.
  ofEdgeId(id: string): { reference: Reference, target: string, paramPath: string } | null {
    const match = EDGE_ID.exec(id);
    if (match === null) {
      return null;
    }
    const kind = this.byMark[match[2]];
    if (kind === undefined) {
      return null;
    }
    return { reference: new Reference(match[1], kind, match[3]), target: match[4], paramPath: match[5] };
  }

  // The kinds a half-typed section still fits.
  couldBe(typedSection: string): ReferenceKind[] {
    return this.kinds.filter(kind => {
      return kind.couldBe(typedSection);
    });
  }

  // Everything one node offers, across all kinds.
  keysOf(node: string, workflowNode: WorkflowNode, catalog: NodeCatalogEntry[]): Reference[] {
    const found: Reference[] = [];
    for (const kind of this.kinds) {
      for (const key of kind.namesOf(workflowNode, catalog)) {
        found.push(new Reference(node, kind, key));
      }
    }
    return found;
  }
}
