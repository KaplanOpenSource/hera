import { flowNodeCenter } from './workflowFlowNodesUtils';

// Cap fit-to-view zoom so a single small node doesn't fill the whole screen.
// Exported for the canvas's own fitView on mount.
export const FIT_MAX_ZOOM = 1;
// How long every viewport move takes.
const MOVE_MS = 300;

// The bits of the ReactFlow instance the viewport needs. Kept narrow so this
// class can be built (and tested) without a canvas.
export interface FlowViewportApi {
  fitView: (options: { duration: number, maxZoom: number }) => void;
  getViewport: () => { x: number, y: number, zoom: number };
  setViewport: (viewport: { x: number, y: number, zoom: number }) => void;
  setCenter: (x: number, y: number, options: { zoom: number, duration: number }) => void;
  getNode: (id: string) => { position: { x: number, y: number }, measured?: { width?: number, height?: number } } | undefined;
}

// Every way the canvas moves itself: fit the whole graph, pan to one node, or
// rescale after the panel's height changed. One place for the zoom cap and the
// animation time, so the moves all feel the same.
export class WorkflowViewport {
  private readonly flow: FlowViewportApi;

  constructor(flow: FlowViewportApi) {
    this.flow = flow;
  }

  // Shows the whole graph.
  fitAll(): void {
    this.flow.fitView({ duration: MOVE_MS, maxZoom: FIT_MAX_ZOOM });
  }

  // Pans to one node, keeping the current zoom. Does nothing when the canvas
  // does not know that node (it may not be measured yet).
  centerOnNode(nodeName: string): void {
    const node = this.flow.getNode(nodeName);
    if (!node) {
      return;
    }
    const center = flowNodeCenter(node);
    this.flow.setCenter(center.x, center.y, { zoom: this.flow.getViewport().zoom, duration: MOVE_MS });
  }

  // Scales the view by `ratio` (the panel grew or shrank in height), so the
  // graph keeps filling the same share of it.
  scaleZoom(ratio: number): void {
    const { x, y, zoom } = this.flow.getViewport();
    this.flow.setViewport({ x: x * ratio, y: y * ratio, zoom: zoom * ratio });
  }
}
