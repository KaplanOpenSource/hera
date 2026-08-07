import { ContentCopy } from '@mui/icons-material';
import { Box, Chip, IconButton, Tooltip } from '@mui/material';
import { pushInfo } from '../../../io/snackbar';
import { LogLineKind } from './classifyLog';
import { KIND_ORDER, KIND_STYLE } from './logKindStyle';

// Filter buttons (one per kind present, with line counts) plus a copy-all button.
export const LogToolbar = ({
  counts,
  visible,
  onToggle,
  fullText,
}: {
  counts: { [kind in LogLineKind]: number },
  visible: { [kind in LogLineKind]: boolean },
  onToggle: (kind: LogLineKind) => void,
  fullText: string,
}) => {
  const copyAll = async () => {
    await navigator.clipboard.writeText(fullText);
    pushInfo('Log copied to clipboard');
  };

  return (
    <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5, flexWrap: 'wrap', mb: 1 }}>
      {KIND_ORDER.filter((kind) => { return counts[kind] > 0; }).map((kind) => {
        return (
          <Chip
            key={kind}
            size="small"
            label={`${KIND_STYLE[kind].label} ${counts[kind]}`}
            variant={visible[kind] ? 'filled' : 'outlined'}
            onClick={() => { return onToggle(kind); }}
            sx={{
              color: KIND_STYLE[kind].color,
              borderColor: KIND_STYLE[kind].color,
              fontWeight: 600,
              opacity: visible[kind] ? 1 : 0.5,
            }}
          />
        );
      })}
      <Box sx={{ flexGrow: 1 }} />
      <Tooltip title="Copy all">
        <IconButton size="small" onClick={copyAll} aria-label="Copy all">
          <ContentCopy fontSize="small" />
        </IconButton>
      </Tooltip>
    </Box>
  );
};
