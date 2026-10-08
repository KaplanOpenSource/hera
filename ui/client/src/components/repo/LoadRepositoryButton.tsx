import { DriveFolderUpload } from "@mui/icons-material";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { loadRepositoryIntoProject } from "../../io/loadRepositoryIntoProject";

export const LoadRepositoryButton = ({
  repositoryName,
}: {
  repositoryName: string,
}) => {
  const [loading, setLoading] = useState(false);

  const handleClick = async () => {
    setLoading(true);
    await loadRepositoryIntoProject(repositoryName);
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
