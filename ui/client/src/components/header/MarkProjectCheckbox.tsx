import { Checkbox } from "@mui/material";

// The box that marks one project in the manage table.
export const MarkProjectCheckbox = ({
  name,
  marked,
  onToggle,
}: {
  name: string,
  marked: boolean,
  onToggle: () => void,
}) => {
  return (
    <Checkbox
      checked={marked}
      onChange={onToggle}
      slotProps={{ input: { 'aria-label': `Mark ${name}` } }}
    />
  );
}
