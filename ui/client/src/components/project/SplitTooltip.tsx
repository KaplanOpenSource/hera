import { Stack, Typography } from "@mui/material";
import { VALUE_GROUP_REST, VALUE_GROUP_UNDEFINED } from "../../utils/splitTree";
import { ToolkitField } from "../details/toolkit/ToolkitField";

// Turns a split value into readable text (handling the special group markers).
const describeValue = (value: string): string => {
  if (value === VALUE_GROUP_REST) return 'other values';
  if (value === VALUE_GROUP_UNDEFINED) return 'not set';
  if (value === '') return "''";
  return value;
};

// Tooltip content for a type / field split node: what the group is, plus its toolkit.
export const SplitTooltip = ({
  fieldLabel,
  value,
  toolkitName,
}: {
  fieldLabel: string,
  value: string,
  toolkitName: string,
}) => {
  const toolkitLabel = toolkitName === VALUE_GROUP_UNDEFINED ? 'No toolkit' : toolkitName;
  return (
    <Stack spacing={1} sx={{ p: 1 }}>
      <Typography variant="subtitle2">Grouped documents</Typography>
      <ToolkitField label="Field" value={fieldLabel} />
      <ToolkitField label="Value" value={describeValue(value)} />
      <ToolkitField label="Toolkit" value={toolkitLabel} />
    </Stack>
  );
};
