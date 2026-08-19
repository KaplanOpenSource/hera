import { FormControlLabel, Switch, Tooltip } from '@mui/material';

// The single view toggle (replaces the old multi-mode selector). OFF shows the
// document's specialized view (agent/workflow editor, or the formulated field
// tree); ON shows the raw stored document with nothing hidden. Always available.
export const RawViewToggle = ({
  rawView,
  setRawView,
}: {
  rawView: boolean,
  setRawView: (v: boolean) => void,
}) => {
  return (
    <Tooltip title="Raw — show the document exactly as stored, with nothing hidden">
      <span>
        <FormControlLabel
          sx={{ mr: 1 }}
          control={
            <Switch
              size="small"
              checked={rawView}
              onChange={(_e, checked) => setRawView(checked)}
            />
          }
          label="Raw View"
        />
      </span>
    </Tooltip>
  );
};
