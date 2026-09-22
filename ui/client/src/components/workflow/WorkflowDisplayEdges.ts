import { Edge, MarkerType } from '@xyflow/react';
import { nodeInputHandleId, nodeOutputHandleId, WorkflowDataflowEdge } from './workflowDataflow';

// Holds the edges ReactFlow draws and builds them up kind by kind. Start with
// the hovered edge id (the X button shows on that one), add the `requires`
// edges and the dataflow edges, then read them back with all():
//
//   WorkflowDisplayEdges.hovering(hoveredEdge)
//     .withRequires(rfEdges, onRemoveRequire)
//     .withDataflow(deps, color, removeDataflowEdge)
//     .all()
//
// Each step returns a new instance, so nothing is edited in place.
export class WorkflowDisplayEdges {
  private readonly edges: Edge[];
  private readonly hoveredEdge: string | null;

  private constructor(edges: Edge[], hoveredEdge: string | null) {
    this.edges = edges;
    this.hoveredEdge = hoveredEdge;
  }

  // An empty set of edges; hoveredEdge is the one the pointer is on, or null.
  static hovering(hoveredEdge: string | null): WorkflowDisplayEdges {
    return new WorkflowDisplayEdges([], hoveredEdge);
  }

  // `requires` edges, attached to the node-level requires handles by id.
  withRequires(
    edges: Edge[],
    onRemove: (source: string, target: string) => void,
  ): WorkflowDisplayEdges {
    const built = edges.map(edge => ({
      ...edge,
      type: 'requires',
      sourceHandle: nodeOutputHandleId(edge.source),
      targetHandle: nodeInputHandleId(edge.target),
      markerEnd: { type: MarkerType.ArrowClosed },
      data: this.dataFor(edge.id, () => onRemove(edge.source, edge.target)),
    }));
    return new WorkflowDisplayEdges([...this.edges, ...built], this.hoveredEdge);
  }

  // Dataflow edges from parameter values that reference another node's output
  // (e.g. `{C.output.ggg}`), drawn output-handle -> input-handle.
  withDataflow(
    edges: WorkflowDataflowEdge[],
    color: string,
    onRemove: (id: string) => void,
  ): WorkflowDisplayEdges {
    const built = edges.map(edge => ({
      ...edge,
      type: 'dataflow',
      markerEnd: { type: MarkerType.ArrowClosed, color },
      style: { stroke: color },
      animated: true,
      // The input handle sits inside the node, so the line's end runs under the
      // node box; lift it above the nodes so it stays visible.
      zIndex: 1000,
      data: this.dataFor(edge.id, () => onRemove(edge.id)),
    }));
    return new WorkflowDisplayEdges([...this.edges, ...built], this.hoveredEdge);
  }

  all(): Edge[] {
    return this.edges;
  }

  // Hover state plus the remove handler the edge's X button calls.
  private dataFor(id: string, onRemove: () => void): { hovered: boolean, onRemove: () => void } {
    return { hovered: id === this.hoveredEdge, onRemove };
  }
}
