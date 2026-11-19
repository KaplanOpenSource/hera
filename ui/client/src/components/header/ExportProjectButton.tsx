import { Download } from "@mui/icons-material";
import { Stack } from "@mui/material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { projectApi } from "../../api/api";
import { DEFAULT_PROJECT, useProjectStore } from "../../stores/useProjectStore";
import { saveFile } from "../../utils/saveFile";

export const ExportProjectButton = ({
  onExported,
}: {
  onExported?: () => void,
}) => {
  const { currProjectName } = useProjectStore();

  const exportProject = async () => {
    const bytes = await projectApi.export(currProjectName);
    if (!bytes) {
      return;
    }
    saveFile(bytes, `${currProjectName}.zip`);
    onExported?.();
  }

  return (
    <ButtonTooltip
      button
      title='Export project'
      disabled={currProjectName === DEFAULT_PROJECT}
      onClick={exportProject}
    >
      <Stack direction={'row'} alignItems={'center'} spacing={0.5}>
        <Download fontSize="small" />
      </Stack>
    </ButtonTooltip>
  );
}
