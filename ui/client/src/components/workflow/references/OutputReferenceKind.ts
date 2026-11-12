import { WorkflowNode } from '../../../shared/types';
import { NodeCatalogEntry, nodeOutputNames } from '../nodeCatalog';
import { ReferenceKind } from './ReferenceKind';

// A reference to another node's output.
export class OutputReferenceKind extends ReferenceKind {
  readonly section = 'output';
  readonly label = 'Output';
  readonly handleMark = 'out';

  namesOf(node: WorkflowNode, catalog: NodeCatalogEntry[]): string[] {
    return nodeOutputNames(node, catalog);
  }

  // Older workflows spell the section `parameters` or `outputs`.
  sectionMatch(): string {
    return 'parameters?|outputs?';
  }
}
