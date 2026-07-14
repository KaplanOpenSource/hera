import { Clear, Search } from '@mui/icons-material';
import { IconButton, InputAdornment, TextField } from '@mui/material';

export const TreeSearchBar = ({
  value,
  onChange,
}: {
  value: string,
  onChange: (value: string) => void,
}) => (
  <TextField
    size="small"
    fullWidth
    placeholder="Search documents…"
    value={value}
    onChange={(e) => onChange(e.target.value)}
    sx={{ mb: 1 }}
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
    }}
  />
);
