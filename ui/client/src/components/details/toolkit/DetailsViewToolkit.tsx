import { Box, Stack, Typography } from '@mui/material';
import { Handyman } from '@mui/icons-material';
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
  const toolkit = toolkitData ? new ToolkitObj(toolkitData) : undefined;
  const documentCount = project?.documents.filter(d => isNoToolkit ? !d.toolkit : d.toolkit === toolkitName).length ?? 0;

  return (
    <Box sx={{ p: 3 }}>
      {toolkit
        ? <ToolkitDetails toolkit={toolkit} documentCount={documentCount} />
        : (
          <Stack spacing={1}>
            <Stack direction="row" spacing={1} alignItems="center">
              <Handyman color="action" />
              <Typography variant="h6">
                {isNoToolkit ? 'Documents without a toolkit' : toolkitName}
              </Typography>
              {!isNoToolkit && (
                <Typography variant="body2" color="text.secondary">(not registered)</Typography>
              )}
            </Stack>
            {documentCount > 0 && (
              <Typography variant="body2" color="text.secondary">
                {documentCount} document{documentCount !== 1 && 's'}
              </Typography>
            )}
          </Stack>
        )
      }
    </Box>
  );
};
