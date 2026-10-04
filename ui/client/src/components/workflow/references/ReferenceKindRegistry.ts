import { WorkflowNode } from '../../../shared/types';
import { NodeCatalogEntry } from '../nodeCatalog';
import { Reference } from './Reference';
import { ReferenceKind } from './ReferenceKind';

// Any reference in a parameter value: the node name, then everything else. Both
// the section and the key may hold dots, so splitSection tells them apart. The
// key may be a JSONPath into the output, e.g. `output.items[0].name`.
const PARSE = /\{\s*(\w+)\.([\w.[\]*'"@$()?,:-]+?)\s*\}/g;
// A source dot's handle id: <node>:<handleMark>:<key>.
const HANDLE = /^(\w+):(\w+):(.+)$/;
// A dataflow edge id: df:<node>:<handleMark>:<key>-><target>.<paramPath>.
// Both the key and the path may hold dots, so `->` is the only split point.
const EDGE_ID = /^df:(\w+):(\w+):(.+?)->(\w+)\.(.+)$/;

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

  // The kind a reference's text after the node name starts with, and the key
  // that follows its section. The longest matching section wins, so a kind whose
  // section holds dots is not cut short by a shorter one. Null when no kind
  // reads it. The key may be empty, for a section that is still being typed.
  splitSection(rest: string): { kind: ReferenceKind, key: string } | null {
    let best: { kind: ReferenceKind, key: string } | null = null;
    for (const kind of this.kinds) {
      const match = new RegExp(`^(?:${kind.sectionMatch()})\\.(.*)$`).exec(rest);
      if (match !== null && (best === null || match[1].length < best.key.length)) {
        best = { kind, key: match[1] };
      }
    }
    return best;
  }

  // Every reference written in one parameter value. Unknown sections are skipped.
  parseAll(value: string): Reference[] {
    const found: Reference[] = [];
    PARSE.lastIndex = 0;
    let match = PARSE.exec(value);
    while (match !== null) {
      const split = this.splitSection(match[2]);
      if (split !== null && split.key !== '') {
        found.push(new Reference(match[1], split.kind, split.key));
      }
      match = PARSE.exec(value);
    }
    return found;
  }

  // One value's text with every known reference to `oldName` pointing at `newName`.
  renamedNode(value: string, oldName: string, newName: string): string {
    return value.replace(PARSE, (match, node, rest) => {
      const split = this.splitSection(rest);
      if (node !== oldName || split === null || split.key === '') {
        return match;
      }
      return `{${newName}.${rest}}`;
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
