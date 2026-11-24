import { Box, LinearProgress, Typography } from "@mui/material";
import { ExportStep } from "./ExportStep";

// A bar and a short line that show how far a multi-project export is.
export const ExportProgress = ({
  step,
}: {
  step: ExportStep,
}) => {
  const value = (step.done / step.total) * 100;

  return (
    <Box sx={{ px: 3, pb: 1 }} data-testid="manage-projects-progress">
      <Typography variant="body2" color="text.secondary" data-testid="manage-projects-progress-text">
        {`Exporting ${step.done + 1} of ${step.total}: ${step.current}`}
      </Typography>
      <LinearProgress
        variant="determinate"
        value={value}
        data-testid="manage-projects-progress-bar"
      />
    </Box>
  );
}
