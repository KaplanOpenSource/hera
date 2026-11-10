import { Box, Typography } from '@mui/material';
import { shownRun, useWorkflowRunStore, WorkflowRunStatus } from '../../../stores/useWorkflowRunStore';
import { WorkflowOutputView } from './WorkflowOutputView';
import { WorkflowRunSelect } from './WorkflowRunSelect';

// One workflow's run output, as its own dock tab. It reads the run straight from
// the shared run store, so it fills live from the poller and survives the run
// button unmounting. A menu button in the bottom bar picks an older saved run instead.
export const WorkflowOutputPanel = ({
  workflowName,
}: {
  workflowName: string,
}) => {
  const run = useWorkflowRunStore((state) => { return shownRun(state, workflowName); });

  const runSelect = <WorkflowRunSelect workflowName={workflowName} />;

  return (
    <Box sx={{ height: '100%', display: 'flex', flexDirection: 'column', minHeight: 0 }}>
      {run && (
        <WorkflowOutputView
          // A different run starts a fresh view, not the old one's scroll spot.
          key={run.token}
          running={run.status === WorkflowRunStatus.Running}
          chunks={run.chunks}
          error={run.status === WorkflowRunStatus.Error ? run.error : null}
          workflowName={workflowName}
          barEnd={runSelect}
        />
      )}
      {!run && (
        <>
          <Typography sx={{ p: 2, flexGrow: 1 }} color="text.secondary">No run yet.</Typography>
          <Box sx={{ display: 'flex', justifyContent: 'flex-end', borderTop: 1, borderColor: 'divider', px: 1, py: 0.5 }}>
            {runSelect}
          </Box>
        </>
      )}
    </Box>
  );
};
