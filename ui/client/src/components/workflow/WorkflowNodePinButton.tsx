import { PushPin, PushPinOutlined } from '@mui/icons-material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';

// The pin among a workflow node's top-right icons: keeps the node's editor open
// after the pointer leaves it. A pinned node stays open until pinned again.
export const WorkflowNodePinButton = ({
  pinned,
  onToggle,
}: {
  pinned: boolean,
  onToggle: () => void,
}) => {
  let icon = <PushPinOutlined sx={{ fontSize: 14 }} />;
  let title = 'Keep it open';
  if (pinned) {
    icon = <PushPin sx={{ fontSize: 14 }} />;
    title = 'Let it close when the pointer leaves';
  }
  return (
    <ButtonTooltip
      className="nodrag"
      aria-label="keep open"
      title={title}
      onClick={onToggle}
      sx={{ p: '2px', bgcolor: 'background.paper', boxShadow: 1 }}
    >
      {icon}
    </ButtonTooltip>
  );
};
