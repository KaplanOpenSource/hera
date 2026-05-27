import { Chip, Stack, Tooltip, Typography } from '@mui/material';
import { Handyman } from '@mui/icons-material';
import { Toolkit } from '../../../shared/types';
import { ToolkitField } from './ToolkitField';

export const ToolkitDetails = ({
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
        {toolkit.source && (
          <Tooltip title="Where this toolkit is registered from (internal or dynamic)">
            <Chip label={toolkit.source} size="small" variant="outlined" />
          </Tooltip>
        )}
        {toolkit.type && (
          <Tooltip title="Toolkit category (measurements, simulations, etc.)">
            <Chip label={toolkit.type} size="small" variant="outlined" />
          </Tooltip>
        )}
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
