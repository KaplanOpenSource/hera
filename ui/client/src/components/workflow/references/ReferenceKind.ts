import { WorkflowNode } from '../../../shared/types';
import { NodeCatalogEntry } from '../nodeCatalog';

// Escapes a string so it can sit inside a regular expression as literal text.
const escapeForRegExp = (text: string): string => {
  return text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
};

// One kind of reference: everything that differs between the kinds of
// `{node.section.key}` a parameter value can hold.
export abstract class ReferenceKind {
  // The written form of the middle part, e.g. `output`.
  abstract readonly section: string;
  // How the kind is named in a menu.
  abstract readonly label: string;
  // This kind's source dot in a handle id, e.g. `out`. One word: handle and edge
  // ids split on dots and colons.
  abstract readonly handleMark: string;

  // The keys of this kind that `node` offers.
  abstract namesOf(node: WorkflowNode, catalog: NodeCatalogEntry[]): string[];

  // The token written into a parameter value.
  format(node: string, key: string): string {
    return `{${node}.${this.section}.${key}}`;
  }

  // The id of the source dot a line of this kind leaves from.
  handleId(node: string, key: string): string {
    return `${node}:${this.handleMark}:${key}`;
  }

  // Matches just this reference, for removing it from a value.
  clearToken(node: string, key: string): RegExp {
    return new RegExp(`\\{\\s*${escapeForRegExp(node)}\\.(?:${this.sectionMatch()})\\.${escapeForRegExp(key)}\\s*\\}`, 'g');
  }

  // The regex source the parser accepts. Override to also read older spellings.
  sectionMatch(): string {
    return escapeForRegExp(this.section);
  }

  // Whether a half-typed section could still become this kind.
  couldBe(typedSection: string): boolean {
    return this.section.startsWith(typedSection);
  }
}
