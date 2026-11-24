import { Checkbox } from "@mui/material";

// The box that marks one project in the manage table.
export const MarkProjectCheckbox = ({
  name,
  marked,
  onToggle,
}: {
  name: string,
  marked: boolean,
  onToggle: (shift: boolean) => void,
}) => {
  return (
    <Checkbox
      checked={marked}
      onChange={(e) => onToggle((e.nativeEvent as MouseEvent).shiftKey)}
      slotProps={{ input: { 'aria-label': `Mark ${name}` } }}
    />
  );
}
