import { Typography } from '@mui/material';
import { useWorkflowRunStore, WorkflowRunStatus } from '../../../stores/useWorkflowRunStore';
import { WorkflowOutputView } from './WorkflowOutputView';

// One workflow's run output, as its own dock tab. It reads the run straight from
// the shared run store, so it fills live from the poller and survives the run
// button unmounting.
export const WorkflowOutputPanel = ({
  workflowName,
}: {
  workflowName: string,
}) => {
  const run = useWorkflowRunStore(state => state.runs[workflowName]);

  return run
    ? (
      <WorkflowOutputView
        running={run.status === WorkflowRunStatus.Running}
        chunks={run.chunks}
        error={run.status === WorkflowRunStatus.Error ? run.error : null}
        workflowName={workflowName}
      />
    )
    : <Typography sx={{ p: 2 }} color="text.secondary">No run yet.</Typography>;
};
