import { Dialog, DialogContent, DialogTitle, Typography } from "@mui/material";
import { ExportProgress } from "./ExportProgress";
import { ExportStep } from "./ExportStep";

// Covers the manage dialog while projects export, so nothing else can be clicked.
export const ExportProgressDialog = ({
  step,
}: {
  step: ExportStep | null,
}) => {
  return (
    <Dialog open fullWidth maxWidth="sm" data-testid="export-progress-dialog">
      <DialogTitle>Exporting projects</DialogTitle>
      <DialogContent>
        {step && <ExportProgress step={step} />}
        {!step && (
          <Typography variant="body2" color="text.secondary">
            Saving projects.zip...
          </Typography>
        )}
      </DialogContent>
    </Dialog>
  );
}
