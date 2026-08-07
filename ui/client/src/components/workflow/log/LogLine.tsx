import { Box } from '@mui/material';
import { ClassifiedLine, LogLineKind } from './classifyLog';

// Per-kind text color and weight. Colors are theme palette keys so they adapt to
// light/dark mode. Output (the task's own stdout) is emphasized on purpose.
const KIND_STYLE: { [kind in LogLineKind]: { color: string, bold?: boolean } } = {
  [LogLineKind.Debug]: { color: 'text.disabled' },
  [LogLineKind.Info]: { color: 'text.secondary' },
  [LogLineKind.Warning]: { color: 'warning.main' },
  [LogLineKind.Error]: { color: 'error.main', bold: true },
  [LogLineKind.Summary]: { color: 'info.main', bold: true },
  [LogLineKind.Output]: { color: 'success.main', bold: true },
};

export const LogLine = ({
  line,
}: {
  line: ClassifiedLine,
}) => {
  const style = KIND_STYLE[line.kind];
  return (
    <Box sx={{ display: 'flex', gap: 1, whiteSpace: 'pre-wrap', wordBreak: 'break-word' }}>
      <Box
        component="span"
        sx={{ color: 'text.disabled', userSelect: 'none', textAlign: 'right', minWidth: 32, flexShrink: 0 }}
      >
        {line.index + 1}
      </Box>
      <Box component="span" sx={{ color: style.color, fontWeight: style.bold ? 600 : 400 }}>
        {line.text || ' '}
      </Box>
    </Box>
  );
};
