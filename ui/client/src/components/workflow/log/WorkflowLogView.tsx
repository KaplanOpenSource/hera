import { Box } from '@mui/material';
import { useLogFilterStore } from '../../../stores/useLogFilterStore';
import { classifyLog } from './classifyLog';
import { LogLine } from './LogLine';

// Renders raw workflow console output as classified, indexed, color-coded lines,
// filtered by the shared log-level visibility. The filter + copy-all toolbar lives
// in the dialog's action bar.
export const WorkflowLogView = ({
  output,
}: {
  output: string,
}) => {
  const visible = useLogFilterStore((state) => { return state.visible; });

  const lines = classifyLog(output);

  return (
    <Box sx={{ fontFamily: 'monospace', fontSize: 12 }}>
      {lines
        .filter((line) => { return visible[line.kind]; })
        .map((line) => { return <LogLine key={line.index} line={line} />; })}
    </Box>
  );
};
