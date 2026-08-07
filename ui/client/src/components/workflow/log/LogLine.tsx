import { Box } from '@mui/material';
import { ClassifiedLine } from './classifyLog';
import { KIND_STYLE } from './logKindStyle';

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
        {line.text || ' '}
      </Box>
    </Box>
  );
};
