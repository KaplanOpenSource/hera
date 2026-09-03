import { Box, Typography } from '@mui/material';
import { WorkflowChunk } from '../../../io/runWorkflow';
import { useLogFilterStore } from '../../../stores/useLogFilterStore';
import { classifyLog, EVENT_PREFIX, LogLineKind } from './classifyLog';
import { LogLine } from './LogLine';
import { LogToolbar } from './LogToolbar';

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

// Drops the [luigi-event] marker lines from a segment's text: in the grouped view
// the task name is already shown as the section label, so the markers are noise.
const withoutEventLines = (
  text: string,
): string => {
  return text
    .split('\n')
    .filter((line) => { return !line.startsWith(EVENT_PREFIX); })
    .join('\n');
};

const emptyCounts = (): { [kind in LogLineKind]: number } => {
  return Object.fromEntries(Object.values(LogLineKind).map((k) => { return [k, 0]; })) as {
    [kind in LogLineKind]: number
  };
};

// Renders a finished run's output grouped by task: each segment is a section with
// a colored left border and the task name (or "setup" / "between") as its label.
// Event marker lines are removed; the log level filters and copy-all still apply.
export const WorkflowChunkedLog = ({
  chunks,
}: {
  chunks: WorkflowChunk[],
}) => {
  const visible = useLogFilterStore((state) => { return state.visible; });
  const toggle = useLogFilterStore((state) => { return state.toggle; });

  // Classify each segment (minus its event lines); keep only segments with content.
  const sections = chunks
    .map((chunk) => {
      const cleaned = withoutEventLines(chunk.text);
      return { name: chunk.name, lines: classifyLog(cleaned) };
    })
    .filter((section) => { return section.lines.some((line) => { return line.text.trim() !== ''; }); });

  const counts = emptyCounts();
  sections.forEach((section) => {
    section.lines.forEach((line) => { counts[line.kind] += 1; });
  });

  const fullText = chunks.map((chunk) => { return withoutEventLines(chunk.text); }).join('');

  return (
    <Box sx={{ fontFamily: 'monospace', fontSize: 12 }}>
      <LogToolbar counts={counts} visible={visible} onToggle={toggle} fullText={fullText} />
      {sections.map((section, sectionIndex) => {
        const style = styleForName(section.name);
        return (
          <Box
            key={sectionIndex}
            sx={{
              borderLeft: 3,
              borderColor: style.color,
              borderLeftStyle: style.dashed ? 'dashed' : 'solid',
              pl: 1,
              ml: 0.5,
              mb: 1.5,
            }}
          >
            <Typography
              variant="caption"
              sx={{ color: style.color, fontWeight: 700, display: 'block', mb: 0.5 }}
            >
              {style.label}
            </Typography>
            {section.lines
              .filter((line) => { return visible[line.kind]; })
              .map((line) => { return <LogLine key={line.index} line={line} />; })}
          </Box>
        );
      })}
    </Box>
  );
};
