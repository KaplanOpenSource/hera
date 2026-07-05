import { Box, Tooltip } from '@mui/material';
import { NodeParameterSource } from '../../shared/types';

// Dot colour + tooltip text for each place a field's value can come from.
const SOURCE_INFO: { [key in NodeParameterSource]: { color: string, description: string } } = {
  [NodeParameterSource.JsonForm]: { color: 'success.main', description: "Declared in the node's jsonForm" },
  [NodeParameterSource.Python]: { color: 'secondary.main', description: "Inferred from the node's Python (executer)" },
  [NodeParameterSource.Template]: { color: 'warning.main', description: "Found in the node's Jinja template" },
};

// A small coloured dot marking where a field's value comes from; hover for the
// explanation. Renders nothing when the field has no known source.
export const FieldSourceDot = ({
  source,
}: {
  source?: NodeParameterSource,
}) => {
  return source ? (
    <Tooltip title={`${source} — ${SOURCE_INFO[source].description}`}>
      <Box sx={{ width: 8, height: 8, borderRadius: '50%', bgcolor: SOURCE_INFO[source].color, flexShrink: 0 }} />
    </Tooltip>
  ) : null;
};
