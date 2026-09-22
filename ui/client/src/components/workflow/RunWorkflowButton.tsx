import { PlayArrow, Save } from '@mui/icons-material';
import { Box, CircularProgress, Menu, MenuItem, SxProps, Theme } from '@mui/material';
import { MouseEvent, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { startWorkflow } from '../../io/runWorkflow';
import { pushError } from '../../io/snackbar';
import { useViewSettingsStore } from '../../stores/useViewSettingsStore';
import { useWorkflowRunStore, WorkflowRunStatus } from '../../stores/useWorkflowRunStore';
import { ProjectDocument } from '../../shared/types';

// Runs a saved workflow via the server. The run happens in the background: starting
// it returns a token, and the shared WorkflowRunPoller polls that token until the
// run finishes. The run lives in useWorkflowRunStore keyed by workflow name, so
// every button for the same workflow shares it: they all disable and show a spinner
// while it runs, and completion/errors reach every one of them. The output shows in
// its own dock tab, opened by the layout when the run enters the store.
//
// Left click runs the workflow. Right click opens a menu with more options:
// run, run with save (one time), and a toggle to always save before running.
// While the toggle is on the icon shows a small save badge and every left click
// saves the document first.
export const RunWorkflowButton = ({
  projectName,
  workflowName,
  doc,
  isChanged,
  save,
  disabled,
  disabledReason,
  sx,
}: {
  projectName: string,
  workflowName: string,
  // The whole workflow document. Sent to the server so it builds from it, no DB lookup.
  doc: ProjectDocument,
  // True when the open document has unsaved edits.
  isChanged?: boolean,
  // Persists the current document; awaited before running when saving is requested.
  save?: () => Promise<void>,
  disabled?: boolean,
  disabledReason?: string,
  sx?: SxProps<Theme>,
}) => {
  const [menuAnchor, setMenuAnchor] = useState<{ x: number, y: number } | null>(null);
  // True only during the brief start request (before the run enters the store).
  const [starting, setStarting] = useState(false);
  const saveBeforeRun = useViewSettingsStore(state => state.viewSettings.alwaysSaveBeforeRun);
  const setViewSettings = useViewSettingsStore(state => state.setViewSettings);
  const run = useWorkflowRunStore(state => state.runs[workflowName]);
  const startRun = useWorkflowRunStore(state => state.startRun);

  const canSave = Boolean(save);
  const isRunning = starting || run?.status === WorkflowRunStatus.Running;

  const doRun = async (withSave: boolean) => {
    setStarting(true);
    try {
      if (withSave && save) {
        await save();
      }
      // Ensure the sent doc carries the resolved name; the server reads it as
      // desc.workflowName (the doc's own desc may leave it unset, falling back to the doc name).
      const docToRun = { ...doc, desc: { ...doc.desc, workflowName } };
      const result = await startWorkflow({ projectName, doc: docToRun });
      if (result.status === 'busy') {
        // A failed start makes no run, so there is no output tab: the snackbar tells it.
        pushError('run workflow: The server is busy running another workflow. Try again shortly.');
      } else if (result.token) {
        // Hand the run to the store; the poller drives it from here.
        startRun(workflowName, result.token);
      }
    } catch (e: any) {
      pushError(`run workflow: ${e?.message ?? e}`);
    } finally {
      setStarting(false);
    }
  };

  // Only save on click when there is something to save.
  const saveOnClick = saveBeforeRun && canSave && Boolean(isChanged);

  const handleClick = () => {
    return doRun(saveOnClick);
  };

  const openMenu = (e: MouseEvent) => {
    e.preventDefault();
    setMenuAnchor({ x: e.clientX, y: e.clientY });
  };

  const closeMenu = () => {
    return setMenuAnchor(null);
  };

  const runFromMenu = (withSave: boolean) => {
    closeMenu();
    doRun(withSave);
  };

  // Disabled while this workflow is running so both buttons block during a run.
  // Unsaved changes no longer block: the run uses the shown workflow, not the saved one.
  const effectiveDisabled = disabled || isRunning;
  let title = 'Run the workflow as shown (right click for options)';
  if (isRunning) {
    title = 'Workflow is running…';
  } else if (disabled && disabledReason) {
    title = disabledReason;
  }
  let icon = <PlayArrow />;
  if (isRunning) {
    icon = <CircularProgress size={18} />;
  } else if (saveBeforeRun && canSave) {
    icon = (
      <Box sx={{ position: 'relative', display: 'inline-flex' }}>
        <PlayArrow />
        <Save sx={{ position: 'absolute', right: -5, bottom: -5, fontSize: 12 }} />
      </Box>
    );
  }

  return (
    <>
      {/* The context menu lives on the wrapper so right click still opens it
          even while the button is disabled (e.g. during a run). */}
      <Box component="span" onContextMenu={openMenu} sx={{ display: 'inline-flex' }}>
        <ButtonTooltip
          title={title}
          aria-label={title}
          disabled={effectiveDisabled}
          onClick={handleClick}
          sx={sx}
        >
          {icon}
        </ButtonTooltip>
      </Box>
      <Menu
        open={Boolean(menuAnchor)}
        onClose={closeMenu}
        anchorReference="anchorPosition"
        anchorPosition={menuAnchor ? { top: menuAnchor.y, left: menuAnchor.x } : undefined}
      >
        <MenuItem onClick={() => runFromMenu(false)} disabled={isRunning}>
          Run as shown
        </MenuItem>
        {canSave && (
          <MenuItem onClick={() => runFromMenu(true)} disabled={isRunning || !isChanged}>
            Run with save
          </MenuItem>
        )}
        {canSave && !saveBeforeRun && (
          <MenuItem onClick={() => { setViewSettings({ alwaysSaveBeforeRun: true }); closeMenu(); }}>
            Auto save before run
          </MenuItem>
        )}
        {canSave && saveBeforeRun && (
          <MenuItem onClick={() => { setViewSettings({ alwaysSaveBeforeRun: false }); closeMenu(); }}>
            Stop auto save before run
          </MenuItem>
        )}
      </Menu>
    </>
  );
};
