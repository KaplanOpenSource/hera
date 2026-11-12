import { OpenInNew } from '@mui/icons-material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';

// The icon among a workflow node's top-right icons that opens the node in a pop-up.
export const WorkflowNodePopupButton = ({
  onOpen,
}: {
  onOpen: () => void,
}) => {
  return (
    <ButtonTooltip
      className="nodrag"
      aria-label="open in a pop-up"
      title="Open in a pop-up"
      onClick={onOpen}
      sx={{ p: '2px', bgcolor: 'background.paper', boxShadow: 1 }}
    >
      <OpenInNew sx={{ fontSize: 14 }} />
    </ButtonTooltip>
  );
};
