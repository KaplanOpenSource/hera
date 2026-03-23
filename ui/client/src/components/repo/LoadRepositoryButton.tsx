import { DriveFolderUpload } from "@mui/icons-material";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { fetchPython } from "../../io/fetchPython";
import { fetchProjectDetails } from "../../io/FetchProjects";
import { useProjectStore } from "../../stores/useProjectStore";

export const LoadRepositoryButton = ({
  repositoryName,
}: {
  repositoryName: string,
}) => {
  const [loading, setLoading] = useState(false);
  const { currProjectName } = useProjectStore();

  const handleClick = async () => {
    setLoading(true);
    await fetchPython({
      results: [],
      code: `
from hera.utils.data.toolkit import dataToolkit
dataToolkit().loadAllDatasourcesInRepositoryToProject(projectName='${currProjectName}', repositoryName='${repositoryName}', overwrite=False)
`,
    });
    await fetchProjectDetails(currProjectName);
    setLoading(false);
  };

  return (
    <ButtonTooltip
      title={`Load repository "${repositoryName}" into project`}
      onClick={handleClick}
      disabled={loading}
    >
      <DriveFolderUpload />
    </ButtonTooltip>
  )
}
