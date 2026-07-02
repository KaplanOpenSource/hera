import { Box, Chip } from '@mui/material';

// A single node output: a chip with the output name and a port dot on the node's
// right edge — the anchor a connection line will leave from (decorative for now;
// becomes a source Handle when outputs get wired). The row is a fixed height so
// it lines up with an input parameter row (whose tree rows don't margin-collapse
// the way a plain stack would).
export const WorkflowNodeOutputChip = ({
  name,
}: {
  name: string,
}) => {
  return (
    <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'flex-end', gap: 0.75, height: '69px' }}>
      <Chip label={name} size="small" variant="outlined" />
      <Box sx={{ width: 8, height: 8, borderRadius: '50%', bgcolor: 'primary.main', flexShrink: 0 }} />
    </Box>
  );
};
