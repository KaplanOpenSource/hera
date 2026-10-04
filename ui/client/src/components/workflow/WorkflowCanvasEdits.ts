import { Connection, Edge } from '@xyflow/react';
import { WorkflowNode } from '../../shared/types';
import { clearInputReference, parseDataflowConnection, parseDataflowEdgeId, setInputReference } from './workflowDataflow';
import { isValidConnection } from './workflowEdges';
import { nodeWithoutParam, nodeWithReferenceAt } from './workflowNodeEdits';
import { Reference } from './references/Reference';

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
      // A node may not reference itself.
      return connection.source !== connection.target;
    }
    return isValidConnection(connection, this.nodeNames, this.nodes);
  }

  // A dragged line: source handle to input handle writes a reference into the
  // target's parameter; otherwise it's requires.
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
      setInputReference(
        this.nodeOf(connection.target),
        dataflow.paramPath,
        new Reference(connection.source, dataflow.kind, dataflow.outputName),
      ),
    );
  }

  // Clears the reference a dataflow edge stands for from its target value.
  removeDataflowEdge(id: string): void {
    const parsed = parseDataflowEdgeId(id);
    if (parsed) {
      this.onSetNode(parsed.target, clearInputReference(this.nodeOf(parsed.target), parsed.paramPath, parsed.refNode, parsed.key));
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

  // Removes one input value from a node (right-click a field -> delete).
  deleteField(nodeName: string, paramPath: string): void {
    this.onSetNode(nodeName, nodeWithoutParam(this.nodeOf(nodeName), paramPath));
  }

  // Inserts a reference into a field's value at the right-click caret.
  referenceOutput(nodeName: string, paramPath: string, reference: Reference, caret?: number): void {
    this.onSetNode(nodeName, nodeWithReferenceAt(this.nodeOf(nodeName), paramPath, reference, caret));
  }

  private nodeOf(name: string): WorkflowNode {
    return this.nodes[name] ?? {};
  }
}
