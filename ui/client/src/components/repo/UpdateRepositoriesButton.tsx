import { Sync } from "@mui/icons-material";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { fetchPython } from "../../io/fetchPython";
import { fetchProjectDetails } from "../../io/FetchProjects";
import { useProjectStore } from "../../stores/useProjectStore";

export const UpdateRepositoriesButton = () => {
  const [loading, setLoading] = useState(false);
  const { currProjectName } = useProjectStore();

  const handleClick = async () => {
    setLoading(true);
    await fetchPython({
      results: [],
      label: 'update repositories',
      code: `
from types import SimpleNamespace
from hera.utils.data.CLI import update
update(SimpleNamespace(projectName='${currProjectName}', overwrite=False))
`,
    });
    await fetchProjectDetails(currProjectName);
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
