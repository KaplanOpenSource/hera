import { Menu, MenuItem } from '@mui/material';

export enum WorkflowContextMenuKind {
  Node = 'node',
  Edge = 'edge',
}

// Right-click target: a node, or an edge (a requires link), anchored at a point.
export type WorkflowContextMenuTarget =
  | { kind: WorkflowContextMenuKind.Node, name: string, x: number, y: number }
  | { kind: WorkflowContextMenuKind.Edge, source: string, target: string, x: number, y: number };

// The right-click menu for the workflow graph: delete a node or remove a
// requires edge, anchored at the click point.
export const WorkflowContextMenu = ({
  menu,
  onClose,
  onDeleteNode,
  onRemoveRequire,
}: {
  menu: WorkflowContextMenuTarget | null,
  onClose: () => void,
  onDeleteNode: (name: string) => void,
  onRemoveRequire: (source: string, target: string) => void,
}) => {
  return (
    <Menu
      open={menu !== null}
      onClose={onClose}
      anchorReference="anchorPosition"
      anchorPosition={menu ? { top: menu.y, left: menu.x } : undefined}
    >
      {menu?.kind === WorkflowContextMenuKind.Node && (
        <MenuItem onClick={() => { onDeleteNode(menu.name); onClose(); }}>
          Delete node “{menu.name}”
        </MenuItem>
      )}
      {menu?.kind === WorkflowContextMenuKind.Edge && (
        <MenuItem onClick={() => { onRemoveRequire(menu.source, menu.target); onClose(); }}>
          Remove requirement ({menu.source} → {menu.target})
        </MenuItem>
      )}
    </Menu>
  );
};
