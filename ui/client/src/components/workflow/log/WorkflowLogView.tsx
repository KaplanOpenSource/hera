import { Box } from '@mui/material';
import { classifyLog } from './classifyLog';
import { LogLine } from './LogLine';

// Renders raw workflow console output as classified, indexed, color-coded lines.
export const WorkflowLogView = ({
  output,
}: {
  output: string,
}) => {
  const lines = classifyLog(output);
  return (
    <Box sx={{ fontFamily: 'monospace', fontSize: 12 }}>
      {lines.map((line) => { return <LogLine key={line.index} line={line} />; })}
    </Box>
  );
};
