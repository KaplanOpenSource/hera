import { Node } from '@xyflow/react';
import { WorkflowNode } from '../../shared/types';
import { NodeCatalogEntry } from './nodeCatalog';
import { NodeRunStatus, NodeRunStatusMap } from './nodeRunStatus';
import { NodePosition } from './WorkflowLayout';
import { computeLayers } from './workflowGeometry';

// Pure helpers for the canvas's node list: rebuilding it after a structure
// change, moving it after a re-layout, the change keys the effects watch, and
// what to do with the viewport once the new nodes are measured.

type Positions = { [id: string]: NodePosition };

// What the canvas should do once the nodes are measured.
export enum FitKind {
  // Fit the whole graph (initial load, or a bulk change like a template).
  All = 'all',
  // Pan to one freshly added node, keeping the current zoom.
  Node = 'node',
  // Nodes only moved: refit right away, no need to wait for measuring.
  Refit = 'refit',
}

export type PendingFit = {
  kind: FitKind,
  // Set only for FitKind.Node.
  nodeName?: string,
};

// The node list rebuilt for the current workflow: dragged positions and
// drag-resized sizes are carried over, new nodes take their laid-out spot.
export const rebuiltFlowNodes = (
  prev: Node[],
  nodeNames: string[],
  layout: Positions,
): Node[] => {
  const prevPos = new Map(prev.map(node => [node.id, node.position]));
  const prevSize = new Map(prev.map(node => [node.id, { width: node.width, height: node.height }]));
  return nodeNames.map(name => ({
    id: name,
    type: 'workflow',
    position: prevPos.get(name) ?? layout[name],
    ...prevSize.get(name),
    data: {},
  }));
};

// The same nodes moved to their laid-out spots (a column changed).
export const restackedFlowNodes = (prev: Node[], layout: Positions): Node[] => {
  return prev.map(node => ({ ...node, position: layout[node.id] ?? node.position }));
};

// The same nodes with only the y positions the de-overlap step changed. Returns
// the list it was given when nothing moved, so React skips the re-render.
export const deOverlappedFlowNodes = (prev: Node[], fixed: Positions): Node[] => {
  let changed = false;
  const next = prev.map(node => {
    const y = fixed[node.id]?.y;
    if (y === undefined || y === node.position.y) {
      return node;
    }
    changed = true;
    return { ...node, position: { ...node.position, y } };
  });
  return changed ? next : prev;
};

// What to fit after the node names changed. A single added node is panned to
// (keeping zoom) so it isn't lost off-screen; the first load and bulk changes
// fit everything; anything else just refits what is there.
export const pendingFitAfterChange = (prevNames: string[], nodeNames: string[]): PendingFit => {
  if (prevNames.length === 0) {
    return { kind: FitKind.All };
  }
  const added = nodeNames.filter(name => !prevNames.includes(name));
  const removed = prevNames.filter(name => !nodeNames.includes(name));
  if (added.length === 0) {
    return { kind: FitKind.Refit };
  }
  if (added.length === 1 && removed.length === 0) {
    return { kind: FitKind.Node, nodeName: added[0] };
  }
  return { kind: FitKind.All };
};

// The middle of a node, from its position and measured size - where the canvas
// centres when it pans to it.
export const flowNodeCenter = (
  node: { position: NodePosition, measured?: { width?: number, height?: number } },
): NodePosition => {
  return {
    x: node.position.x + (node.measured?.width ?? 0) / 2,
    y: node.position.y + (node.measured?.height ?? 0) / 2,
  };
};

// A signature of the workflow structure (names, types, requires, and dataflow
// links) so the graph only rebuilds when the structure changes - not on drag.
export const flowStructureKey = (
  nodeNames: string[],
  nodes: { [name: string]: WorkflowNode },
  dataflowDeps: { source: string, target: string }[],
): string => {
  return JSON.stringify([
    nodeNames.map(name => [name, nodes[name]?.type, nodes[name]?.requires]),
    dataflowDeps.map(edge => [edge.source, edge.target]),
  ]);
};

// A signature of the column each node sits in, so a re-stack happens only when
// a node actually moves column.
export const flowLayerKey = (
  nodeNames: string[],
  nodes: { [name: string]: WorkflowNode },
  dataflowDeps: { source: string, target: string }[],
): string => {
  return JSON.stringify(computeLayers(nodeNames, nodes, dataflowDeps));
};

// A signature of the measured heights, so the de-overlap step runs only after a
// node's real height changed.
export const flowMeasuredKey = (rfNodes: Node[]): string => {
  return JSON.stringify(rfNodes.map(node => [node.id, Math.round(node.measured?.height ?? 0)]));
};

// What a node on the canvas calls back into, each handler naming the node it
// came from (the node component itself only knows its own name).
export interface FlowNodeHandlers {
  onRename: (name: string, newName: string) => void;
  onChange: (name: string, node: WorkflowNode) => void;
  onDelete: (name: string) => void;
  onFieldContextMenu: (name: string, param: string, x: number, y: number, caret?: number) => void;
  onFieldInlineEdit: (name: string, param: string, value: string, caret: number | null, el: HTMLInputElement) => void;
}

// The nodes ReactFlow draws: the canvas's own list with the current selection,
// workflow data and run status laid over it, plus fresh handlers - built each
// render so a node always calls the latest one, with no stale closures.
export const displayFlowNodes = ({
  rfNodes,
  nodes,
  catalog,
  nodeStatuses,
  selectedNode,
  handlers,
}: {
  rfNodes: Node[],
  nodes: { [name: string]: WorkflowNode },
  catalog: NodeCatalogEntry[],
  nodeStatuses?: NodeRunStatusMap,
  selectedNode?: string,
  handlers: FlowNodeHandlers,
}): Node[] => {
  return rfNodes.map(node => ({
    ...node,
    selected: node.id === selectedNode,
    data: {
      name: node.id,
      node: nodes[node.id] ?? {},
      catalog,
      runStatus: nodeStatuses?.[node.id] ?? NodeRunStatus.Pending,
      onRename: (newName: string) => handlers.onRename(node.id, newName),
      onChange: (updated: WorkflowNode) => handlers.onChange(node.id, updated),
      onDelete: () => handlers.onDelete(node.id),
      onFieldContextMenu: (param: string, x: number, y: number, caret?: number) =>
        handlers.onFieldContextMenu(node.id, param, x, y, caret),
      onFieldInlineEdit: (param: string, value: string, caret: number | null, el: HTMLInputElement) =>
        handlers.onFieldInlineEdit(node.id, param, value, caret, el),
    },
  }));
};
