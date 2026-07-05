import { Autocomplete, Box, Menu, MenuItem, TextField } from '@mui/material';
import { ArrowRight } from '@mui/icons-material';
import { ReactNode, useEffect, useRef, useState } from 'react';

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

// The plain (click-to-run) actions for a target. The Field target also gets a
// "Reference output param" item with fly-out submenus, rendered separately.
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

// Props shared by the two fly-out submenus (node picker, then output picker):
// each is a Menu holding a single autocomplete, anchored to its parent item.
const SUBMENU_PROPS = {
  anchorOrigin: { vertical: 'top', horizontal: 'right' },
  transformOrigin: { vertical: 'top', horizontal: 'left' },
  // The autocomplete keeps its own keyboard focus (typing filters rather than
  // driving menu item navigation) and its dropdown overflows the menu paper.
  disableAutoFocusItem: true,
  slotProps: {
    // Tuck each fly-out close to its parent (small negative left) and step it
    // down a touch so the cascade reads as nested; let the dropdown overflow.
    paper: { sx: { overflow: 'visible', ml: -0.5, mt: 0.5 } },
    // Drop the menu list's own top/bottom padding — the box sets its own.
    list: { sx: { py: 0 } },
  },
} as const;

// Padding around each submenu's autocomplete — tight, so the fly-out hugs it.
const SUBMENU_BOX_SX = { px: 1, py: 0.75, width: 220 } as const;

// The right-click menu for the workflow graph: delete a node (or one of its
// input fields), remove a requires edge, or — on a field — reference another
// node's output through two fly-out submenus (pick a node, then pick its output).
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
  // Called once both submenus are chosen: write a reference to sourceNode's
  // output into (node, param).
  onReferenceOutput: (node: string, param: string, sourceNode: string, output: string) => void,
}) => {
  // The item the node submenu flies out from, and the node box the output
  // submenu flies out from.
  const referenceItemRef = useRef<HTMLLIElement | null>(null);
  const nodeBoxRef = useRef<HTMLDivElement | null>(null);
  // Whether the node submenu is open, and which node was chosen in it (which in
  // turn opens the output submenu). Reset when the menu opens or closes.
  const [nodeSubmenuOpen, setNodeSubmenuOpen] = useState(false);
  const [refNode, setRefNode] = useState<string | null>(null);
  useEffect(() => {
    setNodeSubmenuOpen(false);
    setRefNode(null);
  }, [menu]);

  const field = menu?.kind === WorkflowContextMenuKind.Field ? menu : null;
  const showReference = field !== null && referenceOptions.length > 0;
  const outputsForNode = referenceOptions.find(o => o.node === refNode)?.outputs ?? [];
  const actions = menu ? actionsFor(menu, { onDeleteNode, onDeleteField, onRemoveRequire }) : [];

  // Closes the whole cascade (both submenus and the root menu).
  const closeAll = () => {
    setNodeSubmenuOpen(false);
    setRefNode(null);
    onClose();
  };

  return (
    <>
      <Menu
        open={menu !== null}
        onClose={closeAll}
        anchorReference="anchorPosition"
        anchorPosition={menu ? { top: menu.y, left: menu.x } : undefined}
      >
        {showReference && (
          <MenuItem
            ref={referenceItemRef}
            onClick={() => setNodeSubmenuOpen(true)}
            onMouseEnter={() => setNodeSubmenuOpen(true)}
          >
            <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', width: '100%' }}>
              Reference output param
              <ArrowRight fontSize="small" sx={{ ml: 1, mr: -1 }} />
            </Box>
          </MenuItem>
        )}
        {actions.map((action, i) => (
          <MenuItem
            key={i}
            sx={action.danger ? { color: 'error.main' } : undefined}
            onClick={() => { action.run(); closeAll(); }}
          >
            {action.label}
          </MenuItem>
        ))}
        {menu !== null && (
          <MenuItem onClick={closeAll}>Cancel</MenuItem>
        )}
      </Menu>

      {/* Node submenu: pick the source node; choosing one reveals its outputs. */}
      <Menu
        open={showReference && nodeSubmenuOpen}
        anchorEl={referenceItemRef.current}
        onClose={closeAll}
        {...SUBMENU_PROPS}
      >
        <Box ref={nodeBoxRef} sx={SUBMENU_BOX_SX} onKeyDown={e => e.stopPropagation()}>
          <Autocomplete
            size="small"
            openOnFocus
            disablePortal
            value={refNode}
            options={referenceOptions.map(o => o.node)}
            onChange={(_e, value) => setRefNode(value)}
            renderInput={(params) => <TextField {...params} label="Node" autoFocus />}
          />
        </Box>
      </Menu>

      {/* Output submenu (fly-out of the node submenu): pick the output to insert. */}
      <Menu
        open={showReference && nodeSubmenuOpen && refNode !== null}
        anchorEl={nodeBoxRef.current}
        onClose={closeAll}
        {...SUBMENU_PROPS}
      >
        <Box sx={SUBMENU_BOX_SX} onKeyDown={e => e.stopPropagation()}>
          <Autocomplete
            size="small"
            openOnFocus
            disablePortal
            options={outputsForNode}
            onChange={(_e, value) => {
              if (value !== null && field !== null && refNode !== null) {
                onReferenceOutput(field.node, field.param, refNode, value);
                closeAll();
              }
            }}
            renderInput={(params) => <TextField {...params} label="Output" autoFocus />}
          />
        </Box>
      </Menu>
    </>
  );
};
