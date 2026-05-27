import { Box, Stack, Typography } from '@mui/material';
import { FolderOpen, Handyman } from '@mui/icons-material';
import { useProjectStore } from '../../../stores/useProjectStore';
import { ProjectObj } from '../../../objects/ProjectObj';
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

  const toolkit = isNoToolkit ? undefined : toolkits.find(t => t.toolkit === toolkitName || t.shortName === toolkitName);
  const documentCount = project?.documents.filter(d => isNoToolkit ? !d.toolkit : d.toolkit === toolkitName).length ?? 0;

  return (
    <Box sx={{ p: 3 }}>
      {isNoToolkit
        ? (
          <Stack spacing={1}>
            <Stack direction="row" spacing={1} alignItems="center">
              <FolderOpen color="action" />
              <Typography variant="h6">Documents without a toolkit</Typography>
            </Stack>
            <Typography variant="body2" color="text.secondary">
              {documentCount} document{documentCount !== 1 && 's'} not associated with any toolkit
            </Typography>
          </Stack>
        )
        : toolkit
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
