import { Node, useNodesInitialized, useNodesState, useReactFlow } from '@xyflow/react';
import { RefObject, useEffect, useRef } from 'react';
import { WorkflowNode } from '../../shared/types';
import { CanvasResize, canvasResizeAction } from './canvasResize';
import { NodeFocus } from '../../stores/useWorkflowFocusStore';
import { deOverlappedFlowNodes, FitKind, flowLayerKey, flowMeasuredKey, flowStructureKey, PendingFit, pendingFitAfterChange, rebuiltFlowNodes, restackedFlowNodes } from './workflowFlowNodesUtils';
import { WorkflowLayout } from './WorkflowLayout';
import { WorkflowViewport } from './WorkflowViewport';

// Settle time before refitting, so a splitter drag fits once.
const FIT_SETTLE_MS = 120;

// Where the canvas's nodes sit and where the view looks. Keeps the ReactFlow
// node list in step with the workflow (rebuild on a structure change, re-stack
// when a node changes column, push apart nodes that grew) and moves the viewport
// to match (fit after a change, pan to a new or focused node, refit on resize).
// Returns the node list for ReactFlow and the ref for the element to watch.
export const useWorkflowCanvasNodes = ({
  nodeNames,
  nodes,
  dataflowDeps,
  focus,
}: {
  nodeNames: string[],
  nodes: { [name: string]: WorkflowNode },
  dataflowDeps: { source: string, target: string }[],
  focus?: NodeFocus | { nodeName: string, seq: number },
}): {
  rfNodes: Node[],
  onNodesChange: ReturnType<typeof useNodesState<Node>>[2],
  containerRef: RefObject<HTMLDivElement | null>,
} => {
  const viewport = new WorkflowViewport(useReactFlow());
  const nodesInitialized = useNodesInitialized();
  const containerRef = useRef<HTMLDivElement>(null);
  const prevSizeRef = useRef<{ width: number, height: number } | null>(null);
  // What to do with the viewport once the nodes are measured after a structure
  // change. prevNames detects the change.
  const pendingRef = useRef<PendingFit | null>(null);
  const prevNamesRef = useRef<string[]>([]);
  // useNodesState owns only position and identity; the structure effect rebuilds
  // it (preserving dragged positions) when nodes are added/removed/reordered.
  const [rfNodes, setRfNodes, onNodesChange] = useNodesState<Node>([]);

  const structureKey = flowStructureKey(nodeNames, nodes, dataflowDeps);

  // Rebuild the node list when the workflow's structure changes, keeping the
  // positions and sizes the user set on the nodes that survive.
  useEffect(() => {
    const layout = WorkflowLayout.stacked(nodeNames, nodes, dataflowDeps).positions();
    setRfNodes(prev => rebuiltFlowNodes(prev, nodeNames, layout));
  }, [structureKey]);

  // Re-stack when a node moves column: nodes added or removed, or a dependency
  // (requires or a dataflow reference) changed. Not on every drag or param edit.
  const layerKey = flowLayerKey(nodeNames, nodes, dataflowDeps);
  useEffect(() => {
    const layout = WorkflowLayout.stacked(nodeNames, nodes, dataflowDeps).positions();
    setRfNodes(prev => restackedFlowNodes(prev, layout));
    const pending = pendingFitAfterChange(prevNamesRef.current, nodeNames);
    prevNamesRef.current = nodeNames;
    if (pending.kind === FitKind.Refit) {
      requestAnimationFrame(() => viewport.fitAll());
      return;
    }
    // Wait for the new nodes to be measured; fitting now would use stale sizes.
    pendingRef.current = pending;
  }, [layerKey]);

  // Once nodes are measured after a structure change, fit the whole graph or pan
  // to the newly added node.
  useEffect(() => {
    if (!nodesInitialized || pendingRef.current === null) {
      return;
    }
    const pending = pendingRef.current;
    pendingRef.current = null;
    if (pending.kind === FitKind.All) {
      viewport.fitAll();
    } else if (pending.nodeName) {
      viewport.centerOnNode(pending.nodeName);
    }
  }, [nodesInitialized]);

  // The output tab asked for a node: pan to it, keeping the current zoom.
  useEffect(() => {
    if (focus) {
      viewport.centerOnNode(focus.nodeName);
    }
  }, [focus?.seq]);

  // Once nodes are measured, push down only the ones that overlap within their
  // column (using real measured heights) — so growing a node, e.g. by picking a
  // type with more parameters, shoves the nodes below it instead of overlapping
  // them, while leaving every non-colliding position (including drags) untouched.
  const measuredKey = flowMeasuredKey(rfNodes);
  useEffect(() => {
    const fixed = WorkflowLayout.fromFlowNodes(rfNodes, nodeNames, nodes, dataflowDeps).fixOverlaps().positions();
    setRfNodes(prev => deOverlappedFlowNodes(prev, fixed));
  }, [measuredKey]);

  // A width change refits the graph; height alone just scales the zoom to match.
  useEffect(() => {
    const el = containerRef.current;
    if (!el) {
      return;
    }
    let fitTimer: ReturnType<typeof setTimeout>;
    const observer = new ResizeObserver(entries => {
      const { width, height } = entries[0].contentRect;
      const prev = prevSizeRef.current;
      prevSizeRef.current = { width, height };
      const action = canvasResizeAction(prev, { width, height });
      if (action === CanvasResize.Fit) {
        clearTimeout(fitTimer);
        fitTimer = setTimeout(() => viewport.fitAll(), FIT_SETTLE_MS);
      } else if (action === CanvasResize.ScaleZoom) {
        viewport.scaleZoom(height / prev!.height);
      }
    });
    observer.observe(el);
    return () => {
      clearTimeout(fitTimer);
      observer.disconnect();
    };
  }, []);

  return { rfNodes, onNodesChange, containerRef };
};
