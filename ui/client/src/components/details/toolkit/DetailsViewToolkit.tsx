import { Box, Stack, Typography } from '@mui/material';
import { Handyman } from '@mui/icons-material';
import { useProjectStore } from '../../../stores/useProjectStore';
import { ProjectObj } from '../../../objects/ProjectObj';
import { ToolkitDetails } from './ToolkitDetails';

export const DetailsViewToolkit = ({
  toolkitName,
}: {
  toolkitName: string,
}) => {
  const { toolkits, currProject } = useProjectStore();
  const project = currProject ? new ProjectObj(currProject) : null;

  const toolkit = toolkits.find(t => t.toolkit === toolkitName || t.shortName === toolkitName);
  const documentCount = project?.documents.filter(d => d.toolkit === toolkitName).length ?? 0;

  return (
    <Box sx={{ p: 3 }}>
      {toolkit
        ? <ToolkitDetails toolkit={toolkit} documentCount={documentCount} />
        : (
          <Stack direction="row" spacing={1} alignItems="center">
            <Handyman color="action" />
            <Typography variant="h6">{toolkitName}</Typography>
            <Typography variant="body2" color="text.secondary">(not registered)</Typography>
            {documentCount > 0 && (
              <Typography variant="body2" color="text.secondary">{documentCount} document{documentCount !== 1 && 's'}</Typography>
            )}
          </Stack>
        )
      }
    </Box>
  );
};
