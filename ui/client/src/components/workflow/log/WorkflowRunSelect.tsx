import { Menu as MenuIcon } from '@mui/icons-material';
import { Badge, IconButton, Menu, MenuItem, Tooltip } from '@mui/material';
import { MouseEvent, useEffect, useState } from 'react';
import { workflowApi } from '../../../api/api';
import { SavedRunSummary } from '../../../api/WorkflowApi';
import { ProjectDocument } from '../../../shared/types';
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

  // withGraph shows the run's saved graph on the canvas, read-only.
  const pick = async (id: string, withGraph: boolean) => {
    setAnchor(null);
    if (id === LATEST) {
      showLatestRun(workflowName);
      return;
    }
    const run = await workflowApi.run(projectName, workflowName, id);
    if (run) {
      let doc: ProjectDocument | undefined;
      if (withGraph && run.doc) {
        doc = run.doc;
      }
      showOldRun(workflowName, {
        id: run.id,
        token: run.token,
        startTime: run.startTime,
        status: run.status as string as WorkflowRunStatus,
        error: run.error,
        chunks: run.chunks,
        doc,
      });
    }
  };

  // The latest run is already in the list as its own option.
  const older = saved.filter((run) => { return run.token !== latest?.token; });

  // With no run in this session, the newest saved run stands in for the latest, on today's graph.
  useEffect(() => {
    if (!latest && !oldRun && older.length > 0) {
      pick(older[0].id, false);
    }
  }, [latest, oldRun, older.length]);

  // Most users do not know runs are kept, so the tooltip says what the button holds.
  let tooltip = 'Past runs: none yet. Every finished run is saved here.';
  if (older.length === 1) {
    tooltip = 'Past runs: 1 saved run. Click to see its output.';
  } else if (older.length > 1) {
    tooltip = `Past runs: ${older.length} saved runs. Click to see the output of an earlier run.`;
  }

  let value = '';
  if (oldRun) {
    value = oldRun.id;
  } else if (latest) {
    value = LATEST;
  }

  return (
    <>
      <Tooltip title={tooltip}>
        <IconButton
          size="small"
          aria-label="Past runs"
          onClick={(event: MouseEvent<HTMLElement>) => { setAnchor(event.currentTarget); }}
        >
          {/* The number of past runs; hidden when there are none. */}
          <Badge
            badgeContent={older.length}
            color="primary"
            max={99}
            data-testid="past-runs-badge"
            sx={{ '& .MuiBadge-badge': { fontSize: 9, height: 14, minWidth: 14, px: 0.5 } }}
          >
            <MenuIcon fontSize="small" />
          </Badge>
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
          <MenuItem selected={value === LATEST} onClick={() => { pick(LATEST, false); }}>
            {`Latest run - ${latest.status}`}
          </MenuItem>
        )}
        {older.map((run) => {
          return (
            <MenuItem key={run.id} selected={value === run.id} onClick={() => { pick(run.id, true); }}>
              {`${run.startTime} - ${run.status}`}
            </MenuItem>
          );
        })}
        {!latest && older.length === 0 && <MenuItem disabled>No runs yet</MenuItem>}
      </Menu>
    </>
  );
};
