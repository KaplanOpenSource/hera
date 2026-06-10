import { Box, Slider, Typography } from "@mui/material";

// Slider for the auto-reload interval (0–20s). 0 means off, stored as null.
export const ReloadIntervalSlider = ({
  value,
  setValue,
}: {
  value: number | null,
  setValue: (v: number | null) => void,
}) => (
  <Box sx={{ px: 1 }}>
    <Typography variant="caption" color="text.secondary">
      {value === null
        ? 'Auto-reload off (0)'
        : `Auto-reload every ${value} seconds`}
    </Typography>
    <Slider
      value={value ?? 0}
      min={0}
      max={20}
      step={1}
      marks
      valueLabelDisplay="auto"
      onChange={(_, v) => setValue((v as number) === 0 ? null : (v as number))}
    />
  </Box>
);
