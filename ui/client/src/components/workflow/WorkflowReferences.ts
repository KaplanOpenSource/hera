import { WorkflowNode } from '../../shared/types';
import { NodeCatalogEntry } from './nodeCatalog';
import { ReferenceTokenStage, tokenAtCaret } from './workflowDataflow';
import { Reference } from './references/Reference';
import { knownKinds } from './references/knownKinds';

// Answers "which other node's output can this field point at?" for one workflow.
// Build it once from the workflow and the node catalog, then ask it for the
// right-click menu's options or for the suggestions while a {…} reference is
// being typed. Read-only: it never changes a node.
export class WorkflowReferences {
  private readonly nodeNames: string[];
  private readonly nodes: { [name: string]: WorkflowNode };
  private readonly catalog: NodeCatalogEntry[];

  constructor(
    nodeNames: string[],
    nodes: { [name: string]: WorkflowNode },
    catalog: NodeCatalogEntry[],
  ) {
    this.nodeNames = nodeNames;
    this.nodes = nodes;
    this.catalog = catalog;
  }

  // The other nodes that offer something referenceable, with what they offer -
  // what a field on `nodeName` may point at.
  optionsFor(nodeName: string): { node: string, references: Reference[] }[] {
    return this.nodeNames
      .filter(name => name !== nodeName)
      .map(name => ({ node: name, references: this.referencesOf(name) }))
      .filter(option => option.references.length > 0);
  }

  // The suggestions for the `{…}` token the caret sits in, or null when the
  // caret is not inside a token (so the inline menu should close). Node names
  // before the first dot, then the kind labels the typed section still fits,
  // then that kind's keys on the picked node - filtered by the typed text.
  inlineOptions(nodeName: string, value: string, caret: number | null): string[] | null {
    const token = tokenAtCaret(value, caret ?? value.length);
    if (token === null) {
      return null;
    }
    const others = this.optionsFor(nodeName).map(option => option.node);
    const seed = token.seed.toLowerCase();
    if (token.stage === ReferenceTokenStage.Node) {
      return others.filter(name => name.toLowerCase().includes(seed));
    }
    if (!others.includes(token.nodePart)) {
      return [];
    }
    if (token.stage === ReferenceTokenStage.Section) {
      return knownKinds.couldBe(token.sectionPart).map(kind => kind.label);
    }
    return (token.kind?.namesOf(this.nodes[token.nodePart] ?? {}, this.catalog) ?? [])
      .filter(key => key.toLowerCase().includes(seed));
  }

  // Everything one node offers to be referenced.
  referencesOf(nodeName: string): Reference[] {
    return knownKinds.keysOf(nodeName, this.nodes[nodeName] ?? {}, this.catalog);
  }
}
