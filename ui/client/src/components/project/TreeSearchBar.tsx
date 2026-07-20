import { Clear, Search } from '@mui/icons-material';
import { IconButton, InputAdornment, TextField } from '@mui/material';

export const TreeSearchBar = ({
  value,
  onChange,
  warning,
}: {
  value: string,
  onChange: (value: string) => void,
  warning?: string,
}) => (
  <TextField
    size="small"
    fullWidth
    placeholder="Search documents…  (e.g. type:notebook desc.groupName:flat)"
    value={value}
    onChange={(e) => onChange(e.target.value)}
    sx={{ mb: 1 }}
    helperText={warning}
    slotProps={{
      input: {
        startAdornment: (
          <InputAdornment position="start">
            <Search fontSize="small" />
          </InputAdornment>
        ),
        endAdornment: value && (
          <InputAdornment position="end">
            <IconButton size="small" edge="end" onClick={() => onChange('')}>
              <Clear fontSize="small" />
            </IconButton>
          </InputAdornment>
        ),
      },
      formHelperText: { sx: { color: 'warning.main', ml: 0, mt: 0.5 } },
    }}
  />
);
