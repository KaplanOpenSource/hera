import { Add, AutoAwesome } from '@mui/icons-material';
import { Box, Divider, IconButton, Menu, MenuItem, Tooltip, useTheme } from '@mui/material';
import { Background, Connection, Controls, Edge, MarkerType, Node, Panel, ReactFlow, ReactFlowProvider, useNodesInitialized, useNodesState, useReactFlow } from '@xyflow/react';
import '@xyflow/react/dist/style.css';
import { useEffect, useMemo, useRef, useState } from 'react';
import { WorkflowBlock, WorkflowNode } from '../../shared/types';
import { workflowTemplates } from './workflowTemplates';
import { NodeCatalogEntry, nodeOutputNames } from './nodeCatalog';
import { WorkflowContextMenu, WorkflowContextMenuKind, WorkflowContextMenuTarget } from './WorkflowContextMenu';
import { WorkflowFlowNode } from './WorkflowFlowNode';
import { WorkflowRequiresEdge } from './WorkflowRequiresEdge';
import { buildWorkflowEdges, isValidConnection as isValidConnectionPure } from './workflowEdges';
import { buildDataflowEdges, clearInputReference, dataflowReference, insertReferenceAt, nodeInputHandleId, nodeOutputHandleId, parseDataflowConnection, parseDataflowEdgeId, ReferenceTokenStage, replaceReferenceAt, setInputReference, tokenAtCaret } from './workflowDataflow';
import { WorkflowLayout } from './WorkflowLayout';
import { WorkflowInlineReference } from './WorkflowInlineReference';
import { computeLayers } from './workflowGeometry';

// Defined once (module scope) so ReactFlow doesn't warn about changing types.
const NODE_TYPES = { workflow: WorkflowFlowNode };
// Dataflow edges reuse the same removable-edge component (X button at midpoint).
const EDGE_TYPES = { requires: WorkflowRequiresEdge, dataflow: WorkflowRequiresEdge };
// Cap fit-to-view zoom so a single small node doesn't fill the whole screen.
const FIT_MAX_ZOOM = 1;

interface WorkflowGraphProps {
  catalog: NodeCatalogEntry[];
  nodeNames: string[];
  nodes: { [name: string]: WorkflowNode };
  selectedNode?: string;
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
  // Anchor for the templates menu (null when closed).
  const [templatesAnchor, setTemplatesAnchor] = useState<HTMLElement | null>(null);
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
  const prevHeightRef = useRef<number | null>(null);
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

  // When the canvas height changes, scale the zoom by the same ratio so the same
  // slice of the graph stays framed (anchored at the top-left) instead of
  // revealing more or less of it as the height grows or shrinks.
  useEffect(() => {
    const el = containerRef.current;
    if (!el) {
      return;
    }
    const observer = new ResizeObserver(entries => {
      const height = entries[0].contentRect.height;
      const prev = prevHeightRef.current;
      prevHeightRef.current = height;
      if (prev && height && prev !== height) {
        const ratio = height / prev;
        const { x, y, zoom } = getViewport();
        setViewport({ x: x * ratio, y: y * ratio, zoom: zoom * ratio });
      }
    });
    observer.observe(el);
    return () => observer.disconnect();
  }, [getViewport, setViewport]);

  // Overlay current selection and data (with fresh handlers) each render, so the
  // node always calls the latest rename handler — no stale closures, no ref.
  const displayNodes = rfNodes.map(node => ({
    ...node,
    selected: node.id === selectedNode,
    data: {
      name: node.id,
      node: nodes[node.id] ?? {},
      catalog,
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

  // Overlay edge type, direction arrow, hover state, and a fresh remove handler.
  const displayEdges = rfEdges.map(edge => ({
    ...edge,
    type: 'requires',
    // Attach to the node-level requires handles by id.
    sourceHandle: nodeOutputHandleId(edge.source),
    targetHandle: nodeInputHandleId(edge.target),
    markerEnd: { type: MarkerType.ArrowClosed },
    data: {
      hovered: edge.id === hoveredEdge,
      onRemove: () => onRemoveRequire(edge.source, edge.target),
    },
  }));

  // Clears the reference a dataflow edge represents from its target parameter —
  // used both by the edge's X button and by deleting the line.
  const removeDataflowEdge = (id: string) => {
    const ref = parseDataflowEdgeId(id);
    if (ref) {
      onSetNode(ref.target, clearInputReference(nodes[ref.target] ?? {}, ref.param, ref.refNode, ref.key));
    }
  };

  // Dataflow edges from parameter values that reference another node's output
  // (e.g. `{C.output.ggg}`), drawn output-handle → input-handle. Computed each
  // render so edits to parameter values re-derive them. Hovering shows an X that
  // clears the reference, mirroring the requires edges.
  const dataflowEdges: Edge[] = dataflowDeps.map(edge => ({
    ...edge,
    type: 'dataflow',
    markerEnd: { type: MarkerType.ArrowClosed, color: theme.palette.primary.main },
    style: { stroke: theme.palette.primary.main },
    animated: true,
    // The input handle sits inside the node, so the line's end runs under the
    // node box; lift it above the nodes so it stays visible.
    zIndex: 1000,
    data: {
      hovered: edge.id === hoveredEdge,
      onRemove: () => removeDataflowEdge(edge.id),
    },
  }));

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
    const target = nodes[nodeName] ?? {};
    const input_parameters = { ...(target.Execution?.input_parameters ?? {}) };
    delete input_parameters[param];
    onSetNode(nodeName, { ...target, Execution: { ...target.Execution, input_parameters } });
  };

  // While a field's menu is open, the other nodes that produce outputs — the
  // options for its "Reference another node's output…" autocomplete steps.
  const referenceOptions = menu?.kind === WorkflowContextMenuKind.Field
    ? nodeNames
      .filter(name => name !== menu.node)
      .map(name => ({ node: name, outputs: nodeOutputNames(nodes[name] ?? {}, catalog) }))
      .filter(option => option.outputs.length > 0)
    : [];

  // Inserts a {sourceNode.parameters.output} reference into the field's value at
  // the right-click caret (defaulting to the end), leaving the rest of the value
  // intact (the menu closes itself afterward).
  const referenceOutput = (nodeName: string, param: string, sourceNode: string, output: string, caret?: number) => {
    const target = nodes[nodeName] ?? {};
    const params = target.Execution?.input_parameters ?? {};
    const current = typeof params[param] === 'string' ? params[param] : '';
    const next = insertReferenceAt(current, caret ?? current.length, sourceNode, output);
    onSetNode(nodeName, { ...target, Execution: { ...target.Execution, input_parameters: { ...params, [param]: next } } });
  };

  // Other nodes that produce outputs — the nodes a field on `nodeName` can
  // reference (used for both the inline node menu and its output menus).
  const referenceableNodes = (nodeName: string): string[] =>
    nodeNames.filter(name => name !== nodeName).filter(name => nodeOutputNames(nodes[name] ?? {}, catalog).length > 0);

  // The suggestions for the `{…}` token the caret sits in, or null if the caret
  // is not inside a token (so the inline menu should close). Node names before the
  // section dot; the picked node's outputs after it — filtered by the typed text.
  const inlineOptionsFor = (nodeName: string, value: string, caret: number | null): string[] | null => {
    const token = tokenAtCaret(value, caret ?? value.length);
    if (token === null) {
      return null;
    }
    const others = referenceableNodes(nodeName);
    const seed = token.seed.toLowerCase();
    if (token.stage === ReferenceTokenStage.Node) {
      return others.filter(name => name.toLowerCase().includes(seed));
    }
    if (!others.includes(token.nodePart)) {
      return [];
    }
    return nodeOutputNames(nodes[token.nodePart] ?? {}, catalog).filter(output => output.toLowerCase().includes(seed));
  };

  // Typing / caret moves in a field: refresh the inline suggestions, or close them
  // when the caret leaves the token.
  const handleInlineEdit = (nodeName: string, param: string, value: string, caret: number | null, el: HTMLInputElement) => {
    const options = inlineOptionsFor(nodeName, value, caret);
    if (options === null) {
      setInline(null);
      return;
    }
    setInline({ anchorEl: el, node: nodeName, param, options });
  };

  // Writes a new value into the edited field's parameter and queues the caret to
  // land at `caret` once the controlled input re-renders.
  const commitInlineValue = (nodeName: string, param: string, value: string, caret: number, el: HTMLInputElement) => {
    const target = nodes[nodeName] ?? {};
    const params = target.Execution?.input_parameters ?? {};
    onSetNode(nodeName, { ...target, Execution: { ...target.Execution, input_parameters: { ...params, [param]: value } } });
    inlineCaretRef.current = { el, pos: caret };
  };

  // Picks the highlighted inline suggestion. Choosing a node writes the reference
  // scaffold ({node.parameters.}) and switches to picking that node's output;
  // choosing an output completes the {node.parameters.key} token and closes.
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
      const scaffold = `{${option}.parameters.`;
      const next = value.slice(0, token.start) + scaffold + value.slice(token.end);
      commitInlineValue(inline.node, inline.param, next, token.start + scaffold.length, el);
      setInline({ ...inline, options: nodeOutputNames(nodes[option] ?? {}, catalog) });
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
        edges={[...displayEdges, ...dataflowEdges]}
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
        onEdgeContextMenu={(event, edge) => {
          event.preventDefault();
          setMenu({ kind: WorkflowContextMenuKind.Edge, source: edge.source, target: edge.target, x: event.clientX, y: event.clientY });
        }}
      >
        <Panel position="top-right">
          <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1, alignItems: 'flex-end' }}>
            <Tooltip title="Add node">
              <IconButton size="small" onClick={onAddNode} sx={{ bgcolor: 'background.paper', boxShadow: 1, p: 0.5 }}>
                <Add fontSize="small" />
              </IconButton>
            </Tooltip>
            <Tooltip title="Templates">
              <IconButton size="small" onClick={e => setTemplatesAnchor(e.currentTarget)} sx={{ bgcolor: 'background.paper', boxShadow: 1, p: 0.5 }}>
                <AutoAwesome fontSize="small" />
              </IconButton>
            </Tooltip>
          </Box>
          <Menu anchorEl={templatesAnchor} open={!!templatesAnchor} onClose={() => setTemplatesAnchor(null)}>
            {workflowTemplates.map(t => (
              <MenuItem
                key={t.id}
                onClick={() => {
                  onApplyTemplate(t.block);
                  setTemplatesAnchor(null);
                }}
              >
                {t.label}
              </MenuItem>
            ))}
            <Divider />
            <MenuItem onClick={() => setTemplatesAnchor(null)}>Cancel</MenuItem>
          </Menu>
        </Panel>
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
