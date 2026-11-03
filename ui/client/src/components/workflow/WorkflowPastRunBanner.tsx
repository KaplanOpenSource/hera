import { Alert, Box, Button } from '@mui/material';
import { useConfirm } from '../../elements/useConfirm';

// Marks the canvas as showing a past run's saved graph, with a way back to the latest.
export const WorkflowPastRunBanner = ({
  startTime,
  hasUnsavedEdits,
  onRestore,
  onBack,
}: {
  startTime: string,
  // Restore replaces these edits, so it asks first.
  hasUnsavedEdits: boolean,
  onRestore: () => void,
  onBack: () => void,
}) => {
  const { confirmOpen, ConfirmDialog } = useConfirm();

  const restore = async () => {
    if (hasUnsavedEdits) {
      const result = await confirmOpen({
        title: 'Restore this graph?',
        message: 'Your unsaved edits to the workflow will be lost.',
        yesText: 'Restore',
        noText: 'Cancel',
      });
      if (!result.confirmed) {
        return;
      }
    }
    onRestore();
  };

  return (
    <Alert
      severity="info"
      sx={{ py: 0, mb: 1 }}
      action={
        <Box sx={{ display: 'flex', gap: 0.5 }}>
          <Button size="small" color="inherit" onClick={restore}>Restore</Button>
          <Button size="small" color="inherit" onClick={onBack}>Back to latest</Button>
          {ConfirmDialog}
        </Box>
      }
    >
      {`Past run from ${startTime}. This is the graph it ran. Read-only.`}
    </Alert>
  );
};
