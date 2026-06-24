import { Close } from '@mui/icons-material';
import { IconButton } from '@mui/material';

// The little delete "X" shown at the top-right corner of a workflow node.
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
      sx={{ position: 'absolute', top: -12, right: -12, p: '2px', bgcolor: 'background.paper', boxShadow: 1 }}
    >
      <Close sx={{ fontSize: 14 }} />
    </IconButton>
  );
};
