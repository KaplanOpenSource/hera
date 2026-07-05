import { Menu, MenuItem } from '@mui/material';

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
  return (
    <Menu
      open={menu !== null}
      onClose={onClose}
      anchorReference="anchorPosition"
      anchorPosition={menu ? { top: menu.y, left: menu.x } : undefined}
    >
      {menu?.kind === WorkflowContextMenuKind.Field && (
        <MenuItem onClick={() => { onDeleteField(menu.node, menu.param); onClose(); }}>
          Delete input param “{menu.param}”
        </MenuItem>
      )}
      {menu?.kind === WorkflowContextMenuKind.Field && (
        <MenuItem sx={{ color: 'error.main' }} onClick={() => { onDeleteNode(menu.node); onClose(); }}>
          Delete node “{menu.node}”
        </MenuItem>
      )}
      {menu?.kind === WorkflowContextMenuKind.Node && (
        <MenuItem sx={{ color: 'error.main' }} onClick={() => { onDeleteNode(menu.name); onClose(); }}>
          Delete node “{menu.name}”
        </MenuItem>
      )}
      {menu?.kind === WorkflowContextMenuKind.Edge && (
        <MenuItem onClick={() => { onRemoveRequire(menu.source, menu.target); onClose(); }}>
          Remove requirement ({menu.source} → {menu.target})
        </MenuItem>
      )}
      {menu !== null && (
        <MenuItem onClick={onClose}>Cancel</MenuItem>
      )}
    </Menu>
  );
};
