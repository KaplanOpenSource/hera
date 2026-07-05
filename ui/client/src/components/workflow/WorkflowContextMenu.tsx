import { Autocomplete, Box, Menu, MenuItem, TextField } from '@mui/material';
import { ReactNode, useEffect, useState } from 'react';

export enum WorkflowContextMenuKind {
  Node = 'node',
  Field = 'field',
  Edge = 'edge',
}

// Right-click target: a node, one of a node's input fields, or an edge (a
// requires link), anchored at a point.
export type WorkflowContextMenuTarget =
  | { kind: WorkflowContextMenuKind.Node, name: string, x: number, y: number }
  | { kind: WorkflowContextMenuKind.Field, node: string, param: string, x: number, y: number }
  | { kind: WorkflowContextMenuKind.Edge, source: string, target: string, x: number, y: number };

// One node the user can reference, and the outputs it produces.
export interface NodeOutputOption {
  node: string;
  outputs: string[];
}

// One entry in the menu: its label, what it does, and whether it's destructive
// (shown red). `keepOpen` runs the action without closing the menu — used to
// switch the menu into the reference-picking steps in place.
interface MenuAction {
  label: ReactNode,
  run: () => void,
  danger?: boolean,
  keepOpen?: boolean,
}

// The handlers each action can invoke.
interface MenuHandlers {
  onDeleteNode: (name: string) => void,
  onDeleteField: (node: string, param: string) => void,
  onRemoveRequire: (source: string, target: string) => void,
  onStartReference: () => void,
}

// The actions to offer for a given right-click target.
const actionsFor = (menu: WorkflowContextMenuTarget, h: MenuHandlers): MenuAction[] => {
  switch (menu.kind) {
    case WorkflowContextMenuKind.Field:
      return [
        { label: <>Reference another node's output…</>, run: h.onStartReference, keepOpen: true },
        { label: <>Delete input param “{menu.param}”</>, run: () => h.onDeleteField(menu.node, menu.param) },
        { label: <>Delete node “{menu.node}”</>, run: () => h.onDeleteNode(menu.node), danger: true },
      ];
    case WorkflowContextMenuKind.Node:
      return [
        { label: <>Delete node “{menu.name}”</>, run: () => h.onDeleteNode(menu.name), danger: true },
      ];
    case WorkflowContextMenuKind.Edge:
      return [
        { label: <>Remove requirement ({menu.source} → {menu.target})</>, run: () => h.onRemoveRequire(menu.source, menu.target) },
      ];
  }
};

// The right-click menu for the workflow graph: delete a node (or one of its
// input fields), remove a requires edge, or reference another node's output into
// the field — the last via two in-menu autocomplete steps (node, then output).
export const WorkflowContextMenu = ({
  menu,
  referenceOptions,
  onClose,
  onDeleteNode,
  onDeleteField,
  onRemoveRequire,
  onReferenceOutput,
}: {
  menu: WorkflowContextMenuTarget | null,
  // The nodes (and their outputs) offerable when referencing from the current
  // field — only meaningful while a Field menu is open.
  referenceOptions: NodeOutputOption[],
  onClose: () => void,
  onDeleteNode: (name: string) => void,
  onDeleteField: (node: string, param: string) => void,
  onRemoveRequire: (source: string, target: string) => void,
  // Called once both steps are chosen: write a reference to sourceNode's output
  // into (node, param).
  onReferenceOutput: (node: string, param: string, sourceNode: string, output: string) => void,
}) => {
  // The reference flow's step: idle (show actions), pick a node, then pick that
  // node's output. Reset whenever the menu opens on a new target or closes.
  const [refNode, setRefNode] = useState<string | null>(null);
  const [referencing, setReferencing] = useState(false);
  useEffect(() => {
    setReferencing(false);
    setRefNode(null);
  }, [menu]);

  const field = menu?.kind === WorkflowContextMenuKind.Field ? menu : null;
  const actions = menu ? actionsFor(menu, {
    onDeleteNode,
    onDeleteField,
    onRemoveRequire,
    onStartReference: () => setReferencing(true),
  }) : [];
  const outputsForNode = referenceOptions.find(o => o.node === refNode)?.outputs ?? [];

  return (
    <Menu
      open={menu !== null}
      onClose={onClose}
      anchorReference="anchorPosition"
      anchorPosition={menu ? { top: menu.y, left: menu.x } : undefined}
      // The autocomplete step must keep its own keyboard focus; without this the
      // menu grabs it for item type-ahead. Let its dropdown overflow the paper.
      disableAutoFocusItem
      slotProps={{ paper: { sx: { overflow: 'visible' } } }}
    >
      {referencing && field && refNode === null && (
        // Step 1: pick the source node. Stop keydown from reaching the menu so
        // typing filters the autocomplete instead of jumping between items.
        <Box sx={{ px: 2, py: 1, width: 260 }} onKeyDown={e => e.stopPropagation()}>
          <Autocomplete
            size="small"
            openOnFocus
            disablePortal
            options={referenceOptions.map(o => o.node)}
            onChange={(_e, value) => { if (value !== null) { setRefNode(value); } }}
            renderInput={(params) => <TextField {...params} label="Node" autoFocus />}
          />
        </Box>
      )}
      {referencing && field && refNode !== null && (
        // Step 2: pick one of that node's outputs; choosing it inserts and closes.
        <Box sx={{ px: 2, py: 1, width: 260 }} onKeyDown={e => e.stopPropagation()}>
          <Autocomplete
            size="small"
            openOnFocus
            disablePortal
            options={outputsForNode}
            onChange={(_e, value) => {
              if (value !== null) {
                onReferenceOutput(field.node, field.param, refNode, value);
                onClose();
              }
            }}
            renderInput={(params) => <TextField {...params} label={`Output of ${refNode}`} autoFocus />}
          />
        </Box>
      )}
      {!referencing && actions.map((action, i) => (
        <MenuItem
          key={i}
          sx={action.danger ? { color: 'error.main' } : undefined}
          onClick={() => { action.run(); if (!action.keepOpen) { onClose(); } }}
        >
          {action.label}
        </MenuItem>
      ))}
      {menu !== null && !referencing && (
        <MenuItem onClick={onClose}>Cancel</MenuItem>
      )}
    </Menu>
  );
};
