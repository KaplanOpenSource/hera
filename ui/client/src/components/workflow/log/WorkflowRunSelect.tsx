import { Menu as MenuIcon } from '@mui/icons-material';
import { IconButton, Menu, MenuItem, Tooltip } from '@mui/material';
import { MouseEvent, useEffect, useState } from 'react';
import { workflowApi } from '../../../api/api';
import { SavedRunSummary } from '../../../api/WorkflowApi';
import { NO_PROJECT, useProjectStore } from '../../../stores/useProjectStore';
import { useWorkflowRunStore, WorkflowRunStatus } from '../../../stores/useWorkflowRunStore';

// The select value of the run that is going (or just ended) in this session.
const LATEST = '__latest__';

// Picks which run of a workflow the output tab and the canvas show: the latest
// run, or one saved in the run history. A small menu button opens the run list.
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
  const [anchor, setAnchor] = useState<HTMLElement | null>(null);

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
    setAnchor(null);
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
    <>
      <Tooltip title="Pick a run">
        <IconButton
          size="small"
          aria-label="Pick a run"
          onClick={(event: MouseEvent<HTMLElement>) => { setAnchor(event.currentTarget); }}
        >
          <MenuIcon fontSize="small" />
        </IconButton>
      </Tooltip>
      {/* Opens upward from the bar, lined up on its right edge. */}
      <Menu
        anchorEl={anchor}
        open={Boolean(anchor)}
        onClose={() => { setAnchor(null); }}
        anchorOrigin={{ vertical: 'top', horizontal: 'right' }}
        transformOrigin={{ vertical: 'bottom', horizontal: 'right' }}
      >
        {latest && (
          <MenuItem selected={value === LATEST} onClick={() => { pick(LATEST); }}>
            {`Latest run - ${latest.status}`}
          </MenuItem>
        )}
        {older.map((run) => {
          return (
            <MenuItem key={run.id} selected={value === run.id} onClick={() => { pick(run.id); }}>
              {`${run.startTime} - ${run.status}`}
            </MenuItem>
          );
        })}
        {!latest && older.length === 0 && <MenuItem disabled>No runs yet</MenuItem>}
      </Menu>
    </>
  );
};
