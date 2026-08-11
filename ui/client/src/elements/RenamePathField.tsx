import { Box, Tooltip, Typography } from "@mui/material";
import { RenameField } from "./RenameField";

// Like RenameField, but in view mode middle-truncates the path so it shrinks to
// the available space while keeping the last segment visible. Full path in a tooltip.
export const RenamePathField = ({
  value,
  setValue = undefined,
  defaultEditing = false,
}: {
  value: string,
  setValue?: (newVal: string) => void,
  defaultEditing?: boolean,
}) => {
  const splitAt = value.lastIndexOf('/');
  const head = splitAt >= 0 ? value.slice(0, splitAt) : value;
  const tail = splitAt >= 0 ? value.slice(splitAt) : '';

  return (
    <RenameField
      value={value}
      setValue={setValue}
      defaultEditing={defaultEditing}
      valueForView={
        <Tooltip title={value}>
          <Box sx={{ display: 'flex', minWidth: 0, overflow: 'hidden' }}>
            <Typography noWrap sx={{ minWidth: 0 }}>{head}</Typography>
            <Typography sx={{ whiteSpace: 'nowrap', flexShrink: 0 }}>{tail}</Typography>
          </Box>
        </Tooltip>
      }
    />
  )
}
