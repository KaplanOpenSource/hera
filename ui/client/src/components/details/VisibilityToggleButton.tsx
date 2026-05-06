import { Visibility, VisibilityOff } from "@mui/icons-material";
import { Menu, MenuItem } from "@mui/material";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";

export const VisibilityToggleButton = ({
  hasHidden,
  showHidden,
  setShowHidden,
  setHiddenPaths,
  allToggleablePaths,
}: {
  hasHidden: boolean;
  showHidden: boolean;
  setShowHidden: (v: boolean) => void;
  setHiddenPaths: (paths: Set<string>) => void;
  allToggleablePaths: string[];
}) => {
  const [contextMenu, setContextMenu] = useState<{ mouseX: number; mouseY: number } | null>(null);

  return (
    <>
      <ButtonTooltip
        title={showHidden ? 'Show only visible items' : 'Show all items including hidden'}
        onClick={() => setShowHidden(!showHidden)}
        onContextMenu={(e: React.MouseEvent) => {
          e.preventDefault();
          setContextMenu({ mouseX: e.clientX, mouseY: e.clientY });
        }}
      >
        {hasHidden && !showHidden
          ? <VisibilityOff />
          : <Visibility />}
      </ButtonTooltip>
      <Menu
        open={contextMenu !== null}
        onClose={() => setContextMenu(null)}
        anchorReference="anchorPosition"
        anchorPosition={contextMenu ? { top: contextMenu.mouseY, left: contextMenu.mouseX } : undefined}
      >
        <MenuItem onClick={() => { setHiddenPaths(new Set(allToggleablePaths)); setContextMenu(null); }}>
          Hide all
        </MenuItem>
        <MenuItem onClick={() => { setHiddenPaths(new Set()); setContextMenu(null); }}>
          Restore all hidden items
        </MenuItem>
      </Menu>
    </>
  );
};
