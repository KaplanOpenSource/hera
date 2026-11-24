import { Checkbox } from "@mui/material";

// The header box that marks or clears all rows; half-filled when only some are marked.
export const MarkAllCheckbox = ({
  markedCount,
  total,
  onMarkAll,
}: {
  markedCount: number,
  total: number,
  onMarkAll: (mark: boolean) => void,
}) => {
  const allMarked = total > 0 && markedCount === total;
  const someMarked = markedCount > 0 && !allMarked;

  return (
    <Checkbox
      checked={allMarked}
      indeterminate={someMarked}
      onChange={() => onMarkAll(!allMarked)}
      slotProps={{ input: { 'aria-label': 'Mark all' } }}
    />
  );
}
