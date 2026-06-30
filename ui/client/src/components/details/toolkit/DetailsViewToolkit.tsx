import { Box } from '@mui/material';
import { useProjectStore } from '../../../stores/useProjectStore';
import { ProjectObj } from '../../../objects/ProjectObj';
import { ToolkitObj } from '../../../objects/ToolkitObj';
import { ToolkitDetails } from './ToolkitDetails';
import { VALUE_GROUP_UNDEFINED } from '../../../utils/splitTree';

export const DetailsViewToolkit = ({
  toolkitName,
}: {
  toolkitName: string,
}) => {
  const { toolkits, currProject } = useProjectStore();
  const project = currProject ? new ProjectObj(currProject) : null;
  const isNoToolkit = toolkitName === VALUE_GROUP_UNDEFINED;

  const toolkitData = isNoToolkit ? undefined : toolkits.find(t => new ToolkitObj(t).matches(toolkitName));
  const toolkit = toolkitData
    ? new ToolkitObj(toolkitData)
    : ToolkitObj.unregistered(isNoToolkit ? 'Documents without a toolkit' : toolkitName);
  const documentCount = project?.documents.filter(d => toolkit.matches(d.toolkit)).length ?? 0;

  return (
    <Box sx={{ p: 3 }}>
      <ToolkitDetails toolkit={toolkit} documentCount={documentCount} />
    </Box>
  );
};
