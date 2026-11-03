import { Alert, Button } from '@mui/material';

// Marks the canvas as showing a past run's saved graph, with a way back to the latest.
export const WorkflowPastRunBanner = ({
  startTime,
  onBack,
}: {
  startTime: string,
  onBack: () => void,
}) => {
  return (
    <Alert
      severity="info"
      sx={{ py: 0, mb: 1 }}
      action={<Button size="small" color="inherit" onClick={onBack}>Back to latest</Button>}
    >
      {`Past run from ${startTime}. This is the graph it ran. Read-only.`}
    </Alert>
  );
};
