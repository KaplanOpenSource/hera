import { useTheme } from '@mui/material';
import { orange } from '@mui/material/colors';
import { Handle, HandleType, Position } from '@xyflow/react';

export enum RequiresSide {
  DependsOn = 'dependsOn',
  RequiredBy = 'requiredBy',
}

// The big orange dot on a node's left/right edge that you drag to link nodes.
// Left is the "depends on" side, right is the "required by" side.
export const WorkflowRequiresHandle = ({
  id,
  side,
}: {
  id: string,
  side: RequiresSide,
}) => {
  const theme = useTheme();
  let tooltip = 'Depends on: drag to a node that must run before this one';
  let position = Position.Left;
  let handleType: HandleType = 'target';
  if (side === RequiresSide.RequiredBy) {
    tooltip = 'Required by: drag to a node that must run after this one';
    position = Position.Right;
    handleType = 'source';
  }
  // A native title, not MUI's Tooltip: Handle does not forward a ref.
  return (
    <Handle
      type={handleType}
      id={id}
      position={position}
      title={tooltip}
      style={{
        width: 14,
        height: 14,
        background: orange[700],
        border: `2px solid ${theme.palette.background.paper}`,
      }}
    />
  );
};
