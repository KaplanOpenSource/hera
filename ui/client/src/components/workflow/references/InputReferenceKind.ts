import { WorkflowNode } from '../../../shared/types';
import { NodeCatalogEntry, isParameterKey } from '../nodeCatalog';
import { ReferenceKind } from './ReferenceKind';

// A reference to another node's input parameter.
export class InputReferenceKind extends ReferenceKind {
  readonly section = 'Execution.input_parameters';
  readonly label = 'Input';
  // Not `in`: an input row's left dot already uses `<node>:in:<param>`.
  readonly handleMark = 'param';

  namesOf(node: WorkflowNode, _catalog: NodeCatalogEntry[]): string[] {
    return Object.keys(node.Execution?.input_parameters ?? {}).filter(isParameterKey);
  }
}
