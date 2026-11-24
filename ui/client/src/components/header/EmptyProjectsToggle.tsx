import { FormControlLabel, Switch } from "@mui/material";

// Shows only the projects with no documents and no files on disk.
export const EmptyProjectsToggle = ({
  onlyEmpty,
  onChange,
}: {
  onlyEmpty: boolean,
  onChange: (onlyEmpty: boolean) => void,
}) => {
  return (
    <FormControlLabel
      label="Only empty projects"
      control={
        <Switch
          checked={onlyEmpty}
          onChange={(e) => onChange(e.target.checked)}
          slotProps={{ input: { 'aria-label': 'Only empty projects' } }}
        />
      }
    />
  );
}
