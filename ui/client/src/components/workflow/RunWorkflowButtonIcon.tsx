import { PlayArrow, Save, Stop } from '@mui/icons-material';
import { Box, CircularProgress } from '@mui/material';

// The run button's icon: Stop while running, a spinner while starting, else Run.
export const RunWorkflowButtonIcon = ({
  canStop,
  isStarting,
  showSaveBadge,
}: {
  canStop: boolean,
  isStarting: boolean,
  // Shows a small save badge on the Run icon when every run saves first.
  showSaveBadge: boolean,
}) => {
  let icon = <PlayArrow />;
  if (canStop) {
    icon = <Stop />;
  } else if (isStarting) {
    icon = <CircularProgress size={18} />;
  } else if (showSaveBadge) {
    icon = (
      <Box sx={{ position: 'relative', display: 'inline-flex' }}>
        <PlayArrow />
        <Save sx={{ position: 'absolute', right: -5, bottom: -5, fontSize: 12 }} />
      </Box>
    );
  }
  return icon;
};
