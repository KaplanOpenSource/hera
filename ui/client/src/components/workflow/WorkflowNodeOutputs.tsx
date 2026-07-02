import { Box, Typography } from '@mui/material';

// The node's outputs, as a right-aligned column that mirrors the inputs. Each
// row ends in a dot on the node's edge — the anchor a connection line will
// leave from (decorative for now; becomes a source Handle when outputs get
// wired). Rows hide when the inputs tree is collapsed; the title stays.
export const WorkflowNodeOutputs = ({
  outputs,
  expanded,
}: {
  outputs: string[],
  expanded: boolean,
}) => {
  return (
    <Box className="nodrag" sx={{ flexShrink: 0, minWidth: 96, mr: -1 }}>
      <Typography
        sx={{ display: 'flex', alignItems: 'center', justifyContent: 'flex-end', pr: 1, fontSize: '0.875rem', minHeight: '40px' }}
      >
        outputs
      </Typography>
      {expanded && outputs.map(name => (
        <Box
          key={name}
          sx={{ display: 'flex', alignItems: 'center', justifyContent: 'flex-end', gap: 0.75, mt: '7px', mb: '14px', minHeight: '40px' }}
        >
          <Typography sx={{ userSelect: 'text', fontSize: '0.875rem' }}>
            {name}
          </Typography>
          <Box sx={{ width: 8, height: 8, borderRadius: '50%', bgcolor: 'primary.main', flexShrink: 0 }} />
        </Box>
      ))}
    </Box>
  );
};
