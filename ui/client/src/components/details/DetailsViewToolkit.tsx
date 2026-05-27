import { Box, Chip, Stack, Typography } from '@mui/material';
import { Handyman } from '@mui/icons-material';
import { useProjectStore } from '../../stores/useProjectStore';
import { Toolkit } from '../../shared/types';

const ToolkitField = ({
  label,
  value,
}: {
  label: string,
  value: string | undefined,
}) => {
  if (!value) return null;
  return (
    <Stack direction="row" spacing={1} alignItems="baseline">
      <Typography variant="body2" color="text.secondary" sx={{ minWidth: 120 }}>
        {label}
      </Typography>
      <Typography variant="body2" sx={{ fontFamily: 'monospace' }}>
        {value}
      </Typography>
    </Stack>
  );
};

const ToolkitDetails = ({
  toolkit,
  documentCount,
}: {
  toolkit: Toolkit,
  documentCount: number,
}) => {
  return (
    <Stack spacing={2}>
      <Stack direction="row" spacing={1} alignItems="center">
        <Handyman color="action" />
        <Typography variant="h6">{toolkit.toolkit}</Typography>
        {toolkit.source && <Chip label={toolkit.source} size="small" variant="outlined" />}
        {toolkit.type && <Chip label={toolkit.type} size="small" variant="outlined" />}
      </Stack>

      <Stack spacing={1}>
        <ToolkitField label="Class" value={toolkit.cls} />
        <ToolkitField label="Description" value={toolkit.description} />
        <ToolkitField label="Repository" value={toolkit.repositoryName} />
        <ToolkitField label="Version" value={toolkit.version} />
        <ToolkitField label="Documents" value={String(documentCount)} />
      </Stack>
    </Stack>
  );
};

export const DetailsViewToolkit = ({
  toolkitName,
}: {
  toolkitName: string,
}) => {
  const { toolkits } = useProjectStore();
  const project = useProjectStore.getState().getProject();

  const toolkit = toolkits.find(t => t.toolkit === toolkitName);
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
          </Stack>
        )
      }
    </Box>
  );
};
