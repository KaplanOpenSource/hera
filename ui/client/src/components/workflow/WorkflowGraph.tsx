import { Box, useTheme } from '@mui/material';
import { Background, Connection, Controls, Edge, Node, ReactFlow, ReactFlowProvider, useNodesInitialized, useNodesState, useReactFlow } from '@xyflow/react';
import '@xyflow/react/dist/style.css';
import { ReactNode, useEffect, useMemo, useRef, useState } from 'react';
import { WorkflowBlock, WorkflowNode } from '../../shared/types';
import { NodeCatalogEntry } from './nodeCatalog';
import { WorkflowCanvasToolbar } from './WorkflowCanvasToolbar';
import { WorkflowDisplayEdges } from './WorkflowDisplayEdges';
import { nodeWithoutParam, nodeWithParamValue, nodeWithReferenceAt } from './workflowNodeEdits';
import { WorkflowReferences } from './WorkflowReferences';
import { CanvasResize, canvasResizeAction } from './canvasResize';
import { WorkflowContextMenu, WorkflowContextMenuKind, WorkflowContextMenuTarget } from './WorkflowContextMenu';
import { WorkflowFlowNode } from './WorkflowFlowNode';
import { WorkflowRequiresEdge } from './WorkflowRequiresEdge';
import { buildWorkflowEdges, isValidConnection as isValidConnectionPure } from './workflowEdges';
import { buildDataflowEdges, clearInputReference, dataflowReference, parseDataflowConnection, parseDataflowEdgeId, ReferenceTokenStage, replaceReferenceAt, setInputReference, tokenAtCaret } from './workflowDataflow';
import { WorkflowLayout } from './WorkflowLayout';
import { WorkflowInlineReference } from './WorkflowInlineReference';
import { computeLayers } from './workflowGeometry';
import { NodeRunStatus, NodeRunStatusMap } from './nodeRunStatus';

// Defined once (module scope) so ReactFlow doesn't warn about changing types.
const NODE_TYPES = { workflow: WorkflowFlowNode };
// Dataflow edges reuse the same removable-edge component (X button at midpoint).
const EDGE_TYPES = { requires: WorkflowRequiresEdge, dataflow: WorkflowRequiresEdge };
// Cap fit-to-view zoom so a single small node doesn't fill the whole screen.
const FIT_MAX_ZOOM = 1;

// Settle time before refitting, so a splitter drag fits once.
const FIT_SETTLE_MS = 120;

interface WorkflowGraphProps {
  catalog: NodeCatalogEntry[];
  nodeNames: string[];
  nodes: { [name: string]: WorkflowNode };
  selectedNode?: string;
  // How each node did in the last run, keyed by node name. Missing means pending.
  nodeStatuses?: NodeRunStatusMap;
  // A request from the output tab to centre on one node.
  focus?: { nodeName: string, seq: number };
  // The node the pointer is over, or undefined when it left one.
  onHoverNode?: (name: string | undefined) => void;
  // Extra controls rendered in the canvas's top-right corner (e.g. a run button).
  actionButtons?: ReactNode;
  onSelectNode: (name: string | undefined) => void;
  onAddNode: () => void;
  onApplyTemplate: (block: WorkflowBlock) => void;
  onRenameNode: (oldName: string, newName: string) => void;
  onSetNode: (name: string, node: WorkflowNode) => void;
  onAddRequire: (source: string, target: string) => void;
  onRemoveRequire: (source: string, target: string) => void;
  onDeleteNode: (name: string) => void;
}

// Node graph view of a workflow: one node per workflow node, edges from the
// `requires` field. Nodes are draggable, their names editable inline, and the
// on-canvas button adds a node. Clicking a node selects it for editing.
const WorkflowGraphInner = ({
  catalog,
  nodeNames,
  nodes,
  selectedNode,
  nodeStatuses,
  focus,
  onHoverNode,
  actionButtons,
  onSelectNode,
  onAddNode,
  onApplyTemplate,
  onRenameNode,
  onSetNode,
  onAddRequire,
  onRemoveRequire,
  onDeleteNode,
}: WorkflowGraphProps) => {
  const theme = useTheme();
  const [menu, setMenu] = useState<WorkflowContextMenuTarget | null>(null);
  const [hoveredEdge, setHoveredEdge] = useState<string | null>(null);
  // The active inline `{…}` reference autocomplete: which field it hangs under,
  // the node/param being edited, and the current suggestions. Null when idle.
  const [inline, setInline] = useState<{
    anchorEl: HTMLInputElement,
    node: string,
    param: string,
    options: string[],
  } | null>(null);
  // A caret position to restore after an inline pick rewrites the field value
  // (the value is controlled, so we reposition the caret once React re-renders).
  const inlineCaretRef = useRef<{ el: HTMLInputElement, pos: number } | null>(null);
  const { fitView, getViewport, setViewport, getNode, setCenter } = useReactFlow();
  const nodesInitialized = useNodesInitialized();
  const containerRef = useRef<HTMLDivElement>(null);
  const prevSizeRef = useRef<{ width: number, height: number } | null>(null);
  // What to do once nodes are measured after a structure change: 'all' fits the
  // whole graph (initial load / bulk), a node name pans to focus that newly
  // added node while keeping the current zoom. prevNames detects the change.
  const pendingRef = useRef<'all' | string | null>(null);
  const prevNamesRef = useRef<string[]>([]);
  // useNodesState owns only position and identity; the structure effect rebuilds
  // it (preserving dragged positions) when nodes are added/removed/reordered.
  const [rfNodes, setRfNodes, onNodesChange] = useNodesState<Node>([]);

  // Dataflow dependencies (an input referencing another node's output) — used
  // both to draw the lines and to order the columns, alongside `requires`.
  const dataflowDeps = buildDataflowEdges(nodeNames, nodes, catalog);

  // Which node outputs each field may point at (the reference menus ask this).
  const references = new WorkflowReferences(nodeNames, nodes, catalog);

  // A signature of the workflow structure (names, types, requires, and dataflow
  // links) so the graph only rebuilds when the structure changes — not on drag.
  const structureKey = JSON.stringify([
    nodeNames.map(name => [name, nodes[name]?.type, nodes[name]?.requires]),
    dataflowDeps.map(edge => [edge.source, edge.target]),
  ]);

  useEffect(() => {
    const layout = WorkflowLayout.stacked(nodeNames, nodes, dataflowDeps).positions();
    setRfNodes(prev => {
      const prevPos = new Map(prev.map(n => [n.id, n.position]));
      // Carry over any drag-resized dimensions so a structure change (add/remove/
      // retype a node) doesn't reset sizes the user set on the surviving nodes.
      const prevSize = new Map(prev.map(n => [n.id, { width: n.width, height: n.height }]));
      return nodeNames.map(name => ({
        id: name,
        type: 'workflow',
        position: prevPos.get(name) ?? layout[name],
        ...prevSize.get(name),
        data: {},
      }));
    });
  }, [structureKey]);

  // Re-stack whenever nodes are added or removed (and on mount) — but not on
  // every drag or param edit. A single added node is focused (pan to it, keep
  // zoom) so it isn't lost off-screen when fitting all would zoom out too far;
  // initial load / bulk changes fit the whole graph; removes fit what remains.
  // Re-stack when the column assignment changes: nodes added/removed, or a
  // dependency (requires or a dataflow reference) moves a node to another column.
  const layerKey = JSON.stringify(computeLayers(nodeNames, nodes, dataflowDeps));
  useEffect(() => {
    const layout = WorkflowLayout.stacked(nodeNames, nodes, dataflowDeps).positions();
    setRfNodes(prev => prev.map(node => ({ ...node, position: layout[node.id] ?? node.position })));
    const isInitial = prevNamesRef.current.length === 0;
    const added = nodeNames.filter(n => !prevNamesRef.current.includes(n));
    const removed = prevNamesRef.current.filter(n => !nodeNames.includes(n));
    prevNamesRef.current = nodeNames;
    if (isInitial) {
      pendingRef.current = 'all';
    } else if (added.length > 0) {
      // Wait for the new nodes to be measured before fitting. Focus a single
      // freshly-added node (keep zoom); fit everything for bulk/replace changes
      // like applying a template, where fitting now would use unmeasured sizes.
      pendingRef.current = added.length === 1 && removed.length === 0 ? added[0] : 'all';
    } else {
      requestAnimationFrame(() => fitView({ duration: 300, maxZoom: FIT_MAX_ZOOM }));
    }
  }, [layerKey]);

  // Once nodes are measured after a structure change, fit the whole graph or pan
  // to focus the newly added node (keeping the current zoom).
  useEffect(() => {
    if (!nodesInitialized || pendingRef.current === null) {
      return;
    }
    const pending = pendingRef.current;
    pendingRef.current = null;
    if (pending === 'all') {
      fitView({ duration: 300, maxZoom: FIT_MAX_ZOOM });
      return;
    }
    const node = getNode(pending);
    if (node) {
      const x = node.position.x + (node.measured?.width ?? 0) / 2;
      const y = node.position.y + (node.measured?.height ?? 0) / 2;
      setCenter(x, y, { zoom: getViewport().zoom, duration: 300 });
    }
  }, [nodesInitialized]);

  // The output tab asked for a node: pan to it, keeping the current zoom.
  useEffect(() => {
    const node = focus && getNode(focus.nodeName);
    if (!node) {
      return;
    }
    const x = node.position.x + (node.measured?.width ?? 0) / 2;
    const y = node.position.y + (node.measured?.height ?? 0) / 2;
    setCenter(x, y, { zoom: getViewport().zoom, duration: 300 });
  }, [focus?.seq]);

  // Once nodes are measured, push down only the ones that overlap within their
  // column (using real measured heights) — so growing a node, e.g. by picking a
  // type with more parameters, shoves the nodes below it instead of overlapping
  // them, while leaving every non-colliding position (including drags) untouched.
  const measuredKey = JSON.stringify(rfNodes.map(node => [node.id, Math.round(node.measured?.height ?? 0)]));
  useEffect(() => {
    const fixed = WorkflowLayout.fromFlowNodes(rfNodes, nodeNames, nodes, dataflowDeps).fixOverlaps().positions();
    setRfNodes(prev => {
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
    });
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
        fitTimer = setTimeout(() => fitView({ duration: 300, maxZoom: FIT_MAX_ZOOM }), FIT_SETTLE_MS);
      } else if (action === CanvasResize.ScaleZoom) {
        const ratio = height / prev!.height;
        const { x, y, zoom } = getViewport();
        setViewport({ x: x * ratio, y: y * ratio, zoom: zoom * ratio });
      }
    });
    observer.observe(el);
    return () => {
      clearTimeout(fitTimer);
      observer.disconnect();
    };
  }, [getViewport, setViewport, fitView]);

  // Overlay current selection and data (with fresh handlers) each render, so the
  // node always calls the latest rename handler — no stale closures, no ref.
  const displayNodes = rfNodes.map(node => ({
    ...node,
    selected: node.id === selectedNode,
    data: {
      name: node.id,
      node: nodes[node.id] ?? {},
      catalog,
      runStatus: nodeStatuses?.[node.id] ?? NodeRunStatus.Pending,
      onRename: (newName: string) => onRenameNode(node.id, newName),
      onChange: (updated: WorkflowNode) => onSetNode(node.id, updated),
      onDelete: () => onDeleteNode(node.id),
      onFieldContextMenu: (param: string, x: number, y: number, caret?: number) =>
        setMenu({ kind: WorkflowContextMenuKind.Field, node: node.id, param, x, y, caret }),
      onFieldInlineEdit: (param: string, value: string, caret: number | null, el: HTMLInputElement) =>
        handleInlineEdit(node.id, param, value, caret, el),
    },
  }));

  const rfEdges = useMemo<Edge[]>(() => buildWorkflowEdges(nodeNames, nodes), [structureKey]);

  // Clears the reference a dataflow edge represents from its target parameter —
  // used both by the edge's X button and by deleting the line.
  const removeDataflowEdge = (id: string) => {
    const ref = parseDataflowEdgeId(id);
    if (ref) {
      onSetNode(ref.target, clearInputReference(nodes[ref.target] ?? {}, ref.param, ref.refNode, ref.key));
    }
  };

  const displayEdges = WorkflowDisplayEdges.hovering(hoveredEdge)
    .withRequires(rfEdges, onRemoveRequire)
    .withDataflow(dataflowDeps, theme.palette.primary.main, removeDataflowEdge)
    .all();

  const isValidConnection = (connection: Connection | Edge): boolean => {
    // Output→input (dataflow) connections skip the requires cycle check.
    if (parseDataflowConnection(connection.sourceHandle, connection.targetHandle)) {
      return true;
    }
    return isValidConnectionPure(connection, nodeNames, nodes);
  };

  const onConnect = (connection: Connection) => {
    if (!connection.source || !connection.target) {
      return;
    }
    // Dragging an output handle to an input handle writes a dataflow reference
    // ({source.output.name}) into the target's parameter; otherwise it's requires.
    const dataflow = parseDataflowConnection(connection.sourceHandle, connection.targetHandle);
    if (dataflow) {
      onSetNode(connection.target, setInputReference(nodes[connection.target] ?? {}, dataflow.param, connection.source, dataflow.outputName));
      return;
    }
    onAddRequire(connection.source, connection.target);
  };

  // Removes one input parameter from a node (right-click a field → delete).
  const deleteField = (nodeName: string, param: string) => {
    onSetNode(nodeName, nodeWithoutParam(nodes[nodeName] ?? {}, param));
  };

  // While a field's menu is open, the other nodes that produce outputs — the
  // options for its "Reference another node's output…" autocomplete steps.
  let referenceOptions: { node: string, outputs: string[] }[] = [];
  if (menu?.kind === WorkflowContextMenuKind.Field) {
    referenceOptions = references.optionsFor(menu.node);
  }

  // Inserts a {sourceNode.output.name} reference into the field's value at the
  // right-click caret (the menu closes itself afterward).
  const referenceOutput = (nodeName: string, param: string, sourceNode: string, output: string, caret?: number) => {
    onSetNode(nodeName, nodeWithReferenceAt(nodes[nodeName] ?? {}, param, sourceNode, output, caret));
  };

  // Typing / caret moves in a field: refresh the inline suggestions, or close them
  // when the caret leaves the token.
  const handleInlineEdit = (nodeName: string, param: string, value: string, caret: number | null, el: HTMLInputElement) => {
    const options = references.inlineOptions(nodeName, value, caret);
    if (options === null) {
      setInline(null);
      return;
    }
    setInline({ anchorEl: el, node: nodeName, param, options });
  };

  // Writes a new value into the edited field's parameter and queues the caret to
  // land at `caret` once the controlled input re-renders.
  const commitInlineValue = (nodeName: string, param: string, value: string, caret: number, el: HTMLInputElement) => {
    onSetNode(nodeName, nodeWithParamValue(nodes[nodeName] ?? {}, param, value));
    inlineCaretRef.current = { el, pos: caret };
  };

  // Picks the highlighted inline suggestion. Choosing a node writes the reference
  // scaffold ({node.output.}) and switches to picking that node's output;
  // choosing an output completes the {node.output.key} token and closes.
  const pickInline = (option: string) => {
    if (inline === null) {
      return;
    }
    const el = inline.anchorEl;
    const value = el.value;
    const token = tokenAtCaret(value, el.selectionStart ?? value.length);
    if (token === null) {
      setInline(null);
      return;
    }
    if (token.stage === ReferenceTokenStage.Node) {
      const scaffold = `{${option}.output.`;
      const next = value.slice(0, token.start) + scaffold + value.slice(token.end);
      commitInlineValue(inline.node, inline.param, next, token.start + scaffold.length, el);
      setInline({ ...inline, options: references.outputsOf(option) });
    } else {
      const refNode = token.nodePart;
      const next = replaceReferenceAt(value, token.start, token.end, refNode, option);
      commitInlineValue(inline.node, inline.param, next, token.start + dataflowReference(refNode, option).length, el);
      setInline(null);
    }
  };

  // Restore the caret after an inline pick rewrote the (controlled) field value.
  useEffect(() => {
    const pending = inlineCaretRef.current;
    if (pending !== null) {
      inlineCaretRef.current = null;
      pending.el.focus();
      pending.el.setSelectionRange(pending.pos, pending.pos);
    }
  });

  const onEdgesDelete = (deleted: Edge[]) => {
    deleted.forEach(edge => {
      // Deleting a dataflow line clears the reference from its parameter; a
      // requires edge removes the requires link.
      if (parseDataflowEdgeId(edge.id)) {
        removeDataflowEdge(edge.id);
        return;
      }
      onRemoveRequire(edge.source, edge.target);
    });
  };

  return (
    <Box ref={containerRef} sx={{ flex: 1, minHeight: 200, ml: -2, mr: -2, mb: -2, borderTop: '1px solid', borderColor: 'divider' }}>
      <ReactFlow
        colorMode={theme.palette.mode}
        nodes={displayNodes}
        edges={displayEdges}
        nodeTypes={NODE_TYPES}
        edgeTypes={EDGE_TYPES}
        onNodesChange={onNodesChange}
        onConnect={onConnect}
        onEdgesDelete={onEdgesDelete}
        isValidConnection={isValidConnection}
        fitView
        fitViewOptions={{ maxZoom: FIT_MAX_ZOOM }}
        onNodeClick={(_e, node) => onSelectNode(node.id)}
        onPaneClick={() => { onSelectNode(undefined); setInline(null); }}
        onEdgeMouseEnter={(_e, edge) => setHoveredEdge(edge.id)}
        onEdgeMouseLeave={() => setHoveredEdge(null)}
        onNodeContextMenu={(event, node) => {
          event.preventDefault();
          setMenu({ kind: WorkflowContextMenuKind.Node, name: node.id, x: event.clientX, y: event.clientY });
        }}
        onNodeMouseEnter={(_event, node) => { return onHoverNode?.(node.id); }}
        onNodeMouseLeave={() => { return onHoverNode?.(undefined); }}
        onEdgeContextMenu={(event, edge) => {
          event.preventDefault();
          setMenu({ kind: WorkflowContextMenuKind.Edge, source: edge.source, target: edge.target, x: event.clientX, y: event.clientY });
        }}
      >
        <WorkflowCanvasToolbar
          actionButtons={actionButtons}
          onAddNode={onAddNode}
          onApplyTemplate={onApplyTemplate}
        />
        <Background />
        <Controls />
      </ReactFlow>
      <WorkflowContextMenu
        menu={menu}
        referenceOptions={referenceOptions}
        onClose={() => setMenu(null)}
        onDeleteNode={onDeleteNode}
        onDeleteField={deleteField}
        onRemoveRequire={onRemoveRequire}
        onReferenceOutput={referenceOutput}
      />
      <WorkflowInlineReference
        anchorEl={inline?.anchorEl ?? null}
        options={inline?.options ?? []}
        onPick={pickInline}
        onClose={() => setInline(null)}
      />
    </Box>
  );
};

// ReactFlowProvider supplies the store that hooks like useNodesInitialized read,
// so the inner component must live inside it.
export const WorkflowGraph = (props: WorkflowGraphProps) => {
  return (
    <ReactFlowProvider>
      <WorkflowGraphInner {...props} />
    </ReactFlowProvider>
  );
};
