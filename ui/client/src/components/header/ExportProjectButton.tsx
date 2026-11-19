import { Download } from "@mui/icons-material";
import { Stack } from "@mui/material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { projectApi } from "../../api/api";
import { DEFAULT_PROJECT, useProjectStore } from "../../stores/useProjectStore";

// Makes the browser save the bytes as a file on the user's machine.
const saveFile = (bytes: ArrayBuffer, fileName: string) => {
  const url = URL.createObjectURL(new Blob([bytes], { type: 'application/zip' }));
  const link = document.createElement('a');
  link.href = url;
  link.download = fileName;
  link.click();
  URL.revokeObjectURL(url);
}

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
