import { Box, Tooltip } from '@mui/material';
import { NodeParameterSource } from '../../shared/types';

// Dot colour + tooltip text for each place a field's value can come from.
const SOURCE_INFO: { [key in NodeParameterSource]: { color: string, description: string } } = {
  [NodeParameterSource.JsonForm]: { color: 'success.main', description: "Declared in the node's jsonForm" },
  [NodeParameterSource.Python]: { color: 'secondary.main', description: "Inferred from the node's Python (executer)" },
  [NodeParameterSource.Template]: { color: 'warning.main', description: "Found in the node's Jinja template" },
};

// The dot's colour and tooltip for a given source, the grey "unknown" dot when
// there is no source but the caller opted into showing one, or null (no dot).
const dotInfo = (
  source?: NodeParameterSource,
  showUnknown?: boolean,
): { color: string, title: string } | null => {
  if (source) {
    return { color: SOURCE_INFO[source].color, title: `${source} — ${SOURCE_INFO[source].description}` };
  }
  if (showUnknown) {
    return { color: 'action.disabled', title: "No declared source — drag another node's output here to connect" };
  }
  return null;
};

// A small coloured dot marking where a field's value comes from; hover for the
// explanation. Renders nothing when the field has no known source — unless
// `showUnknown` is set, in which case an unsourced field gets a grey dot (used
// on workflow nodes so it can still be a drag target for another node's output).
export const FieldSourceDot = ({
  source,
  showUnknown = false,
}: {
  source?: NodeParameterSource,
  showUnknown?: boolean,
}) => {
  const info = dotInfo(source, showUnknown);
  return info ? (
    <Tooltip title={info.title}>
      <Box sx={{ width: 8, height: 8, borderRadius: '50%', bgcolor: info.color, flexShrink: 0 }} />
    </Tooltip>
  ) : null;
};
