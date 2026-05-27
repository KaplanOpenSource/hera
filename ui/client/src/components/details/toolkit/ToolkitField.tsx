import { Stack, Tooltip, Typography } from '@mui/material';

const FIELD_TOOLTIPS: Record<string, string> = {
  'Class': 'Fully-qualified Python class path of the toolkit implementation',
  'Description': 'Short description of what the toolkit does',
  'Repository': 'Repository JSON file where this toolkit is defined',
  'Version': 'Toolkit version',
  'Documents': 'Number of documents in this project belonging to this toolkit',
};

export const ToolkitField = ({
  label,
  value,
}: {
  label: string,
  value: string | undefined,
}) => {
  if (!value) return null;
  return (
    <Stack direction="row" spacing={1} alignItems="baseline">
      <Tooltip title={FIELD_TOOLTIPS[label] ?? ''} placement="left">
        <Typography variant="body2" color="text.secondary" sx={{ minWidth: 120 }}>
          {label}
        </Typography>
      </Tooltip>
      <Typography variant="body2" sx={{ fontFamily: 'monospace' }}>
        {value}
      </Typography>
    </Stack>
  );
};
