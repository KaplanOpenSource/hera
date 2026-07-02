import { Box, Typography } from '@mui/material';
import { WorkflowNodeOutputChip } from './WorkflowNodeOutputChip';

// The node's outputs, as a right-aligned column that mirrors the inputs: one
// chip per output, each on a row the height of an input row so they line up.
// Rows hide when the inputs tree is collapsed; the title stays.
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
        // Center in a box the height of the inputs' title row (its small icon
        // buttons) so the two titles line up.
        sx={{ display: 'flex', alignItems: 'center', justifyContent: 'flex-end', pr: 1, fontSize: '0.875rem', minHeight: '38px' }}
      >
        outputs
      </Typography>
      {expanded && outputs.map(name => (
        <WorkflowNodeOutputChip key={name} name={name} />
      ))}
    </Box>
  );
};
