import { Button } from "@mui/material";
import { useNavigate } from "react-router-dom";
import { useProjectStore } from "../../stores/useProjectStore";

// Makes this project the current one, the same way the project dropdown does.
export const SwitchToProjectButton = ({
  name,
  onSwitched,
}: {
  name: string,
  onSwitched: () => void,
}) => {
  const { currProjectName } = useProjectStore();
  const navigate = useNavigate();
  const isCurrent = name === currProjectName;

  const switchTo = () => {
    navigate('/' + encodeURIComponent(name));
    onSwitched();
  }

  let label = 'Switch';
  if (isCurrent) {
    label = 'Current';
  }

  return (
    <Button
      size="small"
      disabled={isCurrent}
      onClick={switchTo}
      aria-label={`Switch to ${name}`}
      data-testid="manage-projects-switch"
    >
      {label}
    </Button>
  );
}
