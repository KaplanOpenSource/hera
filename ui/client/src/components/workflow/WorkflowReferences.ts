import { WorkflowNode } from '../../shared/types';
import { NodeCatalogEntry, nodeOutputNames } from './nodeCatalog';
import { ReferenceTokenStage, tokenAtCaret } from './workflowDataflow';

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

  // The other nodes that produce outputs, with those outputs - what a field on
  // `nodeName` may reference.
  optionsFor(nodeName: string): { node: string, outputs: string[] }[] {
    return this.nodeNames
      .filter(name => name !== nodeName)
      .map(name => ({ node: name, outputs: this.outputsOf(name) }))
      .filter(option => option.outputs.length > 0);
  }

  // The suggestions for the `{…}` token the caret sits in, or null when the
  // caret is not inside a token (so the inline menu should close). Node names
  // before the section dot; the picked node's outputs after it - filtered by
  // the typed text.
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
    return this.outputsOf(token.nodePart).filter(output => output.toLowerCase().includes(seed));
  }

  // The outputs one node produces.
  outputsOf(nodeName: string): string[] {
    return nodeOutputNames(this.nodes[nodeName] ?? {}, this.catalog);
  }
}
