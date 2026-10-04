import { WorkflowNode } from '../../shared/types';
import { computeLayers, estimateHeight, estimatedWidth, H_GAP, V_GAP, X_GAP } from './workflowGeometry';

// A node placed on the canvas: the column (dependency layer) it belongs to, its
// position, and its size.
export interface PlacedNode {
  id: string;
  layer: number;
  x: number;
  y: number;
  height: number;
  width: number;
}

export interface NodePosition {
  x: number;
  y: number;
}

// Holds the positions and sizes of the workflow's nodes and can de-overlap them.
// Build it from a workflow (estimated heights), from the canvas's live nodes
// (real measured heights), or from explicit placed nodes; then read back
// positions(). Stacking a fresh layout and fixing overlaps after a node grew are
// the same operation — push each column's nodes down until they no longer
// collide — so both go through fixOverlaps().
export class WorkflowLayout {
  private readonly placed: PlacedNode[];

  private constructor(placed: PlacedNode[]) {
    this.placed = placed;
  }

  static fromPlaced(placed: PlacedNode[]): WorkflowLayout {
    return new WorkflowLayout(placed);
  }

  // Fresh layout from a workflow: every node at the top of its column with an
  // estimated height, then stacked so they do not overlap.
  static stacked(
    nodeNames: string[],
    nodes: { [name: string]: WorkflowNode },
    extraDeps: { source: string, target: string }[] = [],
  ): WorkflowLayout {
    const layers = computeLayers(nodeNames, nodes, extraDeps);
    const placed = nodeNames.map(id => {
      const layer = layers[id] ?? 0;
      return { id, layer, x: layer * X_GAP, y: 0, height: estimateHeight(nodes[id] ?? {}), width: estimatedWidth() };
    });
    return new WorkflowLayout(placed).fixOverlaps();
  }

  // Layout from the canvas's current nodes, keeping their positions and using
  // their real measured heights (falling back to an estimate before a node is
  // measured). Pair with fixOverlaps() to resolve overlaps after a size change.
  static fromFlowNodes(
    flowNodes: { id: string, position: NodePosition, measured?: { height?: number, width?: number } }[],
    nodeNames: string[],
    nodes: { [name: string]: WorkflowNode },
    extraDeps: { source: string, target: string }[] = [],
  ): WorkflowLayout {
    const layers = computeLayers(nodeNames, nodes, extraDeps);
    const placed = flowNodes.map(node => ({
      id: node.id,
      layer: layers[node.id] ?? 0,
      x: node.position.x,
      y: node.position.y,
      height: node.measured?.height ?? estimateHeight(nodes[node.id] ?? {}),
      width: node.measured?.width ?? estimatedWidth(),
    }));
    return new WorkflowLayout(placed);
  }

  // Push down only the nodes that overlap within their column — never up or
  // sideways, never the topmost node — so a node that grew shoves the ones below
  // it just enough to clear, while every non-colliding position (including
  // deliberate drags) is preserved. Mutates this layout in place and returns it.
  fixOverlaps(vGap: number = V_GAP): this {
    Object.values(this.columns()).forEach(column => {
      let prevBottom = -Infinity;
      column
        .slice()
        .sort((a, b) => a.y - b.y)
        .forEach(node => {
          node.y = Math.max(node.y, prevBottom);
          prevBottom = node.y + node.height + vGap;
        });
    });
    return this;
  }

  // Pull the nodes together, using the sizes this layout was built with: each
  // column sits one gap right of the widest node in the column before it, and
  // within a column the topmost node stays where it is while every node below
  // sits one gap under the one above. Built from the canvas's measured nodes,
  // this is what keeps the canvas as tight as the nodes really are - and what
  // makes a node that grew, e.g. one pinned open, push the others aside.
  // Mutates this layout in place and returns it.
  compact(vGap: number = V_GAP, hGap: number = H_GAP): this {
    this.spreadColumns(hGap);
    Object.values(this.columns()).forEach(column => {
      let prevBottom: number | null = null;
      column
        .slice()
        .sort((a, b) => a.y - b.y)
        .forEach(node => {
          if (prevBottom !== null) {
            node.y = prevBottom;
          }
          prevBottom = node.y + node.height + vGap;
        });
    });
    return this;
  }

  // Put each column one gap right of the widest node in the column before it.
  // The leftmost column stays where it is.
  private spreadColumns(hGap: number): void {
    const columns = this.columns();
    let x: number | null = null;
    Object.keys(columns)
      .map(Number)
      .sort((a, b) => a - b)
      .forEach(layer => {
        const column = columns[layer];
        if (x === null) {
          x = Math.min(...column.map(node => node.x));
        }
        const left = x;
        column.forEach(node => {
          node.x = left;
        });
        x = left + Math.max(...column.map(node => node.width)) + hGap;
      });
  }

  // The nodes of each column (dependency layer), keyed by layer.
  private columns(): { [layer: number]: PlacedNode[] } {
    const columns: { [layer: number]: PlacedNode[] } = {};
    this.placed.forEach(node => {
      (columns[node.layer] ??= []).push(node);
    });
    return columns;
  }

  // The {x, y} of every node, keyed by id.
  positions(): { [id: string]: NodePosition } {
    const result: { [id: string]: NodePosition } = {};
    this.placed.forEach(node => {
      result[node.id] = { x: node.x, y: node.y };
    });
    return result;
  }
}
