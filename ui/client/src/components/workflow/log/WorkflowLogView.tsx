import { Box } from '@mui/material';
import { useState } from 'react';
import { classifyLog, LogLineKind } from './classifyLog';
import { LogLine } from './LogLine';
import { LogToolbar } from './LogToolbar';

// All kinds start visible.
const allVisible = (): { [kind in LogLineKind]: boolean } => {
  return Object.fromEntries(Object.values(LogLineKind).map((k) => { return [k, true]; })) as {
    [kind in LogLineKind]: boolean
  };
};

const emptyCounts = (): { [kind in LogLineKind]: number } => {
  return Object.fromEntries(Object.values(LogLineKind).map((k) => { return [k, 0]; })) as {
    [kind in LogLineKind]: number
  };
};

// Renders raw workflow console output as classified, indexed, color-coded lines,
// with per-kind show/hide filters (showing counts) and a copy-all button.
export const WorkflowLogView = ({
  output,
}: {
  output: string,
}) => {
  const [visible, setVisible] = useState(allVisible);

  const lines = classifyLog(output);

  const counts = emptyCounts();
  lines.forEach((line) => { counts[line.kind] += 1; });

  const toggle = (kind: LogLineKind) => {
    setVisible((prev) => { return { ...prev, [kind]: !prev[kind] }; });
  };

  return (
    <Box sx={{ fontFamily: 'monospace', fontSize: 12 }}>
      <LogToolbar counts={counts} visible={visible} onToggle={toggle} fullText={output} />
      {lines
        .filter((line) => { return visible[line.kind]; })
        .map((line) => { return <LogLine key={line.index} line={line} />; })}
    </Box>
  );
};
