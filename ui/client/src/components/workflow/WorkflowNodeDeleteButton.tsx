import { Close } from '@mui/icons-material';
import { IconButton } from '@mui/material';

// The delete "X" among a workflow node's top-right icons. The node's icon slot
// places it; this only draws it.
export const WorkflowNodeDeleteButton = ({
  onDelete,
}: {
  onDelete: () => void,
}) => {
  return (
    <IconButton
      className="nodrag"
      size="small"
      onClick={(e) => { e.stopPropagation(); onDelete(); }}
      sx={{ p: '2px', bgcolor: 'background.paper', boxShadow: 1 }}
    >
      <Close sx={{ fontSize: 14 }} />
    </IconButton>
  );
};
