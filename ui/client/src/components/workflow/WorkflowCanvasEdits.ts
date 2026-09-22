import { Connection, Edge } from '@xyflow/react';
import { WorkflowNode } from '../../shared/types';
import { clearInputReference, parseDataflowConnection, parseDataflowEdgeId, setInputReference } from './workflowDataflow';
import { isValidConnection } from './workflowEdges';
import { nodeWithoutParam, nodeWithReferenceAt } from './workflowNodeEdits';

// What a gesture on the canvas does to the workflow: drawing a line, deleting
// one, and the field commands on the right-click menu. Each one works out the
// new node and hands it to the callbacks the editor passed in - nothing here
// keeps state of its own.
export class WorkflowCanvasEdits {
  private readonly nodeNames: string[];
  private readonly nodes: { [name: string]: WorkflowNode };
  private readonly onSetNode: (name: string, node: WorkflowNode) => void;
  private readonly onAddRequire: (source: string, target: string) => void;
  private readonly onRemoveRequire: (source: string, target: string) => void;

  constructor(
    nodeNames: string[],
    nodes: { [name: string]: WorkflowNode },
    callbacks: {
      onSetNode: (name: string, node: WorkflowNode) => void,
      onAddRequire: (source: string, target: string) => void,
      onRemoveRequire: (source: string, target: string) => void,
    },
  ) {
    this.nodeNames = nodeNames;
    this.nodes = nodes;
    this.onSetNode = callbacks.onSetNode;
    this.onAddRequire = callbacks.onAddRequire;
    this.onRemoveRequire = callbacks.onRemoveRequire;
  }

  // Output→input (dataflow) connections skip the requires cycle check.
  canConnect(connection: Connection | Edge): boolean {
    if (parseDataflowConnection(connection.sourceHandle, connection.targetHandle)) {
      return true;
    }
    return isValidConnection(connection, this.nodeNames, this.nodes);
  }

  // A dragged line: output handle to input handle writes a dataflow reference
  // ({source.output.name}) into the target's parameter; otherwise it's requires.
  connect(connection: Connection): void {
    if (!connection.source || !connection.target) {
      return;
    }
    const dataflow = parseDataflowConnection(connection.sourceHandle, connection.targetHandle);
    if (!dataflow) {
      this.onAddRequire(connection.source, connection.target);
      return;
    }
    this.onSetNode(
      connection.target,
      setInputReference(this.nodeOf(connection.target), dataflow.param, connection.source, dataflow.outputName),
    );
  }

  // Clears the reference a dataflow edge stands for from its target parameter.
  removeDataflowEdge(id: string): void {
    const ref = parseDataflowEdgeId(id);
    if (ref) {
      this.onSetNode(ref.target, clearInputReference(this.nodeOf(ref.target), ref.param, ref.refNode, ref.key));
    }
  }

  // Deleting lines: a dataflow line clears its reference, a requires line drops
  // the requires link.
  removeEdges(edges: Edge[]): void {
    edges.forEach(edge => {
      if (parseDataflowEdgeId(edge.id)) {
        this.removeDataflowEdge(edge.id);
        return;
      }
      this.onRemoveRequire(edge.source, edge.target);
    });
  }

  // Removes one input parameter from a node (right-click a field -> delete).
  deleteField(nodeName: string, param: string): void {
    this.onSetNode(nodeName, nodeWithoutParam(this.nodeOf(nodeName), param));
  }

  // Inserts a {sourceNode.output.name} reference into a field's value at the
  // right-click caret.
  referenceOutput(nodeName: string, param: string, sourceNode: string, output: string, caret?: number): void {
    this.onSetNode(nodeName, nodeWithReferenceAt(this.nodeOf(nodeName), param, sourceNode, output, caret));
  }

  private nodeOf(name: string): WorkflowNode {
    return this.nodes[name] ?? {};
  }
}
