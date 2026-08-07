import { Box } from '@mui/material';
import { useLogFilterStore } from '../../../stores/useLogFilterStore';
import { classifyLog, LogLineKind } from './classifyLog';
import { LogLine } from './LogLine';
import { LogToolbar } from './LogToolbar';

const emptyCounts = (): { [kind in LogLineKind]: number } => {
  return Object.fromEntries(Object.values(LogLineKind).map((k) => { return [k, 0]; })) as {
    [kind in LogLineKind]: number
  };
};

// Renders raw workflow console output as classified, indexed, color-coded lines,
// with per-kind show/hide filters (showing counts) and a copy-all button. The
// filter choices are persisted (useLogFilterStore) so they survive across runs.
export const WorkflowLogView = ({
  output,
}: {
  output: string,
}) => {
  const visible = useLogFilterStore((state) => { return state.visible; });
  const toggle = useLogFilterStore((state) => { return state.toggle; });

  const lines = classifyLog(output);

  const counts = emptyCounts();
  lines.forEach((line) => { counts[line.kind] += 1; });

  return (
    <Box sx={{ fontFamily: 'monospace', fontSize: 12 }}>
      <LogToolbar counts={counts} visible={visible} onToggle={toggle} fullText={output} />
      {lines
        .filter((line) => { return visible[line.kind]; })
        .map((line) => { return <LogLine key={line.index} line={line} />; })}
    </Box>
  );
};
