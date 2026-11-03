import { MenuItem, TextField } from '@mui/material';
import { useEffect, useState } from 'react';
import { workflowApi } from '../../../api/api';
import { SavedRunSummary } from '../../../api/WorkflowApi';
import { NO_PROJECT, useProjectStore } from '../../../stores/useProjectStore';
import { useWorkflowRunStore, WorkflowRunStatus } from '../../../stores/useWorkflowRunStore';

// The select value of the run that is going (or just ended) in this session.
const LATEST = '__latest__';

// Picks which run of a workflow the output tab and the canvas show: the latest
// run, or one saved in the run history.
export const WorkflowRunSelect = ({
  workflowName,
}: {
  workflowName: string,
}) => {
  const projectName = useProjectStore((state) => { return state.currProjectName; });
  const latest = useWorkflowRunStore((state) => { return state.runs[workflowName]; });
  const oldRun = useWorkflowRunStore((state) => { return state.oldRuns[workflowName]; });
  const showOldRun = useWorkflowRunStore((state) => { return state.showOldRun; });
  const showLatestRun = useWorkflowRunStore((state) => { return state.showLatestRun; });
  const [saved, setSaved] = useState<SavedRunSummary[]>([]);

  // Read the list again when a run starts or ends, so the run before it shows up.
  useEffect(() => {
    if (projectName === NO_PROJECT) {
      return;
    }
    let cancelled = false;
    const load = async () => {
      const runs = await workflowApi.runs(projectName, workflowName);
      if (!cancelled) {
        setSaved(runs ?? []);
      }
    };
    load();
    return () => {
      cancelled = true;
    };
  }, [projectName, workflowName, latest?.token, latest?.status]);

  const pick = async (id: string) => {
    if (id === LATEST) {
      showLatestRun(workflowName);
      return;
    }
    const run = await workflowApi.run(projectName, workflowName, id);
    if (run) {
      showOldRun(workflowName, {
        id: run.id,
        token: run.token,
        status: run.status as string as WorkflowRunStatus,
        error: run.error,
        chunks: run.chunks,
      });
    }
  };

  // The latest run is already in the list as its own option.
  const older = saved.filter((run) => { return run.token !== latest?.token; });

  // With no run in this session, show the newest saved run instead of nothing.
  useEffect(() => {
    if (!latest && !oldRun && older.length > 0) {
      pick(older[0].id);
    }
  }, [latest, oldRun, older.length]);

  let value = '';
  if (oldRun) {
    value = oldRun.id;
  } else if (latest) {
    value = LATEST;
  }

  return (
    <TextField
      select
      size="small"
      label="Run"
      value={value}
      onChange={(event) => { pick(event.target.value); }}
      sx={{ minWidth: 260 }}
    >
      {latest && <MenuItem value={LATEST}>{`Latest run - ${latest.status}`}</MenuItem>}
      {older.map((run) => {
        return <MenuItem key={run.id} value={run.id}>{`${run.startTime} - ${run.status}`}</MenuItem>;
      })}
    </TextField>
  );
};
