import { Menu, MenuItem } from '@mui/material';
import { ReactNode } from 'react';

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

// One entry in the menu: its label, what it does, and whether it's destructive
// (shown red).
interface MenuAction {
  label: ReactNode,
  run: () => void,
  danger?: boolean,
}

// The handlers each action can invoke.
interface MenuHandlers {
  onDeleteNode: (name: string) => void,
  onDeleteField: (node: string, param: string) => void,
  onRemoveRequire: (source: string, target: string) => void,
}

// The actions to offer for a given right-click target.
const actionsFor = (menu: WorkflowContextMenuTarget, h: MenuHandlers): MenuAction[] => {
  switch (menu.kind) {
    case WorkflowContextMenuKind.Field:
      return [
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
// input fields) or remove a requires edge, anchored at the click point.
export const WorkflowContextMenu = ({
  menu,
  onClose,
  onDeleteNode,
  onDeleteField,
  onRemoveRequire,
}: {
  menu: WorkflowContextMenuTarget | null,
  onClose: () => void,
  onDeleteNode: (name: string) => void,
  onDeleteField: (node: string, param: string) => void,
  onRemoveRequire: (source: string, target: string) => void,
}) => {
  const actions = menu ? actionsFor(menu, { onDeleteNode, onDeleteField, onRemoveRequire }) : [];
  return (
    <Menu
      open={menu !== null}
      onClose={onClose}
      anchorReference="anchorPosition"
      anchorPosition={menu ? { top: menu.y, left: menu.x } : undefined}
    >
      {actions.map((action, i) => (
        <MenuItem
          key={i}
          sx={action.danger ? { color: 'error.main' } : undefined}
          onClick={() => { action.run(); onClose(); }}
        >
          {action.label}
        </MenuItem>
      ))}
      {menu !== null && (
        <MenuItem onClick={onClose}>Cancel</MenuItem>
      )}
    </Menu>
  );
};
