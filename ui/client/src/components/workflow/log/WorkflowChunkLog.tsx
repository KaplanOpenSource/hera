import { Paper, Typography } from '@mui/material';
import { WorkflowChunk } from '../../../io/runWorkflow';
import { useLogFilterStore } from '../../../stores/useLogFilterStore';
import { classifyLog } from './classifyLog';
import { LogLine } from './LogLine';
import { withoutEventLines } from './logMetrics';

// Segment names the server uses for output that is not tied to a task.
const PREAMBLE = '__preamble__';
const BETWEEN = '__between__';

type SegmentStyle = {
  label: string,
  color: string,
  dashed: boolean,
};

const styleForName = (
  name: string,
): SegmentStyle => {
  if (name === PREAMBLE) {
    return { label: 'setup', color: 'text.disabled', dashed: true };
  }
  if (name === BETWEEN) {
    return { label: 'between', color: 'text.disabled', dashed: true };
  }
  return { label: name, color: 'primary.main', dashed: false };
};

// Renders one output chunk (a single task's segment) as a self-contained boxed card:
// the task name (or "setup" / "between") as a header, then its classified lines.
// Reads the shared log-level filter so it stays in sync with the toolbar. Kept on
// its own so it can later show a single task's log on the canvas. Renders nothing
// when the chunk has no visible content.
export const WorkflowChunkLog = ({
  chunk,
}: {
  chunk: WorkflowChunk,
}) => {
  const visible = useLogFilterStore((state) => { return state.visible; });

  const style = styleForName(chunk.name);
  const lines = classifyLog(withoutEventLines(chunk.text));
  const hasContent = lines.some((line) => { return line.text.trim() !== ''; });

  return hasContent && (
    <Paper
      variant="outlined"
      sx={{
        fontFamily: 'monospace',
        fontSize: 12,
        p: 1,
        mb: 1.5,
        borderColor: style.color,
        borderStyle: style.dashed ? 'dashed' : 'solid',
      }}
    >
      <Typography
        variant="caption"
        sx={{ color: style.color, fontWeight: 700, display: 'block', mb: 0.5 }}
      >
        {style.label}
      </Typography>
      {lines
        .filter((line) => { return visible[line.kind]; })
        .map((line) => { return <LogLine key={line.index} line={line} />; })}
    </Paper>
  );
};
