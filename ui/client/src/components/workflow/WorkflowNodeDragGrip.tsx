import { DragIndicator } from '@mui/icons-material';
import { Box } from '@mui/material';

// The grab bar at the top middle of a node. The open editor swallows drags, so
// this is the one spot that always moves the node.
export const WorkflowNodeDragGrip = () => {
  return (
    <Box
      title="Drag to move this node"
      sx={{
        position: 'absolute',
        top: -11,
        left: '50%',
        transform: 'translateX(-50%)',
        display: 'flex',
        alignItems: 'center',
        px: 0.5,
        borderRadius: 1,
        bgcolor: 'background.paper',
        border: '1px solid',
        borderColor: 'divider',
        color: 'text.secondary',
        cursor: 'grab',
        zIndex: 1002,
      }}
    >
      <DragIndicator sx={{ fontSize: 14, transform: 'rotate(90deg)' }} />
    </Box>
  );
};
