import { Autocomplete, TextField } from '@mui/material';

// Multi-select for a node's `requires` field, constrained to the given node
// names. Normalizes the stored value (a single name or a list) to an array.
export const WorkflowRequiresSelect = ({
  requires,
  options,
  onChange,
}: {
  requires?: string | string[],
  options: string[],
  onChange: (values: string[]) => void,
}) => {
  const value = requires === undefined
    ? []
    : (Array.isArray(requires) ? requires : [requires]);

  return (
    <Autocomplete
      multiple
      size="small"
      sx={{ minWidth: 220 }}
      options={options}
      value={value}
      onChange={(_e, values) => onChange(values)}
      renderInput={(params) => <TextField {...params} label="requires" />}
    />
  );
};
