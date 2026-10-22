import { Sync } from "@mui/icons-material";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { loadAllRepositoriesIntoProject } from "../../io/loadRepositoryIntoProject";

export const UpdateRepositoriesButton = () => {
  const [loading, setLoading] = useState(false);

  const handleClick = async () => {
    setLoading(true);
    // The answer brings the project along, so there is nothing to reload after it.
    await loadAllRepositoriesIntoProject();
    setLoading(false);
  };

  return (
    <ButtonTooltip
      title="Update repositories into project"
      onClick={handleClick}
      disabled={loading}
    >
      <Sync />
    </ButtonTooltip>
  )
}
