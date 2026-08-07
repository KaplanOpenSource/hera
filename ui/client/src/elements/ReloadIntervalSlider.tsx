import { Box, Slider, Stack, Typography } from "@mui/material";
import { AutoReloadToggle } from "../components/header/AutoReloadToggle";

// Slider for the auto-reload interval (0–20s). 0 means off, stored as null.
// The on/off toggle sits on the right of the title (same control as the toolbar).
export const ReloadIntervalSlider = ({
  value,
  setValue,
}: {
  value: number | null,
  setValue: (v: number | null) => void,
}) => (
  <Box>
    <Typography>
      Auto-reload
      <Typography component="span" variant="caption" color="text.secondary" sx={{ ml: 0.5 }}>
        {value === null ? 'off (0)' : `every ${value} seconds`}
      </Typography>
    </Typography>
    <Stack direction="row" spacing={1} alignItems="center" sx={{ mb: 2 }}>
      <AutoReloadToggle />
      <Slider
        value={value ?? 0}
        min={0}
        max={20}
        step={1}
        marks
        valueLabelDisplay="auto"
        onChange={(_, v) => setValue((v as number) === 0 ? null : (v as number))}
      />
    </Stack>
  </Box>
);
