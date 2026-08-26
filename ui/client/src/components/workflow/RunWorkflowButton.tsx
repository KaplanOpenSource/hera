import { PlayArrow, Save } from '@mui/icons-material';
import { Box, Menu, MenuItem, SxProps, Theme } from '@mui/material';
import { MouseEvent, useEffect, useRef, useState } from 'react';
import { SnackbarKey } from 'notistack';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { startWorkflow, pollWorkflow } from '../../io/runWorkflow';
import { dismiss, pushError, pushInfo, pushRunning } from '../../io/snackbar';
import { useViewSettingsStore } from '../../stores/useViewSettingsStore';
import { WorkflowOutputDialog } from './log/WorkflowOutputDialog';

// How often to poll a running workflow for its status.
const POLL_MS = 500;

// Runs a saved workflow via the server. The run happens in the background: starting
// it returns a token, and the button polls that token every POLL_MS until the run
// is done. The output dialog opens right away and shows a spinner until the
// captured console output (or an error) arrives.
//
// Left click runs the workflow. Right click opens a menu with more options:
// run, run with save (one time), and a toggle to always save before running.
// While the toggle is on the icon shows a small save badge and every left click
// saves the document first.
export const RunWorkflowButton = ({
  projectName,
  workflowName,
  isChanged,
  save,
  disabled,
  disabledReason,
  sx,
}: {
  projectName: string,
  workflowName: string,
  // True when the open document has unsaved edits.
  isChanged?: boolean,
  // Persists the current document; awaited before running when saving is requested.
  save?: () => Promise<void>,
  disabled?: boolean,
  disabledReason?: string,
  sx?: SxProps<Theme>,
}) => {
  const [open, setOpen] = useState(false);
  const [running, setRunning] = useState(false);
  const [output, setOutput] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  // The token of the run currently being polled; null when nothing is polling.
  const [token, setToken] = useState<string | null>(null);
  const [menuAnchor, setMenuAnchor] = useState<{ x: number, y: number } | null>(null);
  const saveBeforeRun = useViewSettingsStore(state => state.viewSettings.alwaysSaveBeforeRun);
  const setViewSettings = useViewSettingsStore(state => state.setViewSettings);
  // Key of the "run workflow" running-snackbar, so we can dismiss it when the run ends.
  const runningKeyRef = useRef<SnackbarKey | null>(null);

  const canSave = Boolean(save);

  // Clears the running state and dismisses the running snackbar. Setting token to
  // null also stops the poll effect.
  const finishRun = () => {
    setRunning(false);
    setToken(null);
    if (runningKeyRef.current) {
      dismiss(runningKeyRef.current);
      runningKeyRef.current = null;
    }
  };

  // Poll the server for the run's status until it is done or failed. Modeled on
  // ServerReadyGate: a setTimeout loop with a cancelled guard, cleaned up on
  // unmount (or when the token changes).
  useEffect(() => {
    if (!token) {
      return;
    }
    let cancelled = false;
    let timer: ReturnType<typeof setTimeout>;

    const poll = async () => {
      let result;
      try {
        result = await pollWorkflow(token);
      } catch (e: any) {
        if (cancelled) {
          return;
        }
        setError(e?.message ?? String(e));
        pushError(`run workflow: ${e?.message ?? e}`);
        finishRun();
        return;
      }
      if (cancelled) {
        return;
      }
      if (result.status === 'running') {
        timer = setTimeout(poll, POLL_MS);
      } else if (result.status === 'done') {
        setOutput(result.output ?? '');
        pushInfo(`Workflow "${workflowName}" finished`);
        finishRun();
      } else if (result.status === 'error') {
        const message = result.error || 'Workflow failed';
        setError(message);
        pushError(`run workflow: ${message}`);
        finishRun();
      } else {
        // not_found: the server restarted or a newer run overwrote the slot.
        const message = 'The run was lost (the server may have restarted).';
        setError(message);
        pushError(`run workflow: ${message}`);
        finishRun();
      }
    };
    poll();

    return () => {
      cancelled = true;
      clearTimeout(timer);
    };
  }, [token]);

  const doRun = async (withSave: boolean) => {
    setOpen(true);
    setRunning(true);
    setOutput(null);
    setError(null);
    setToken(null);
    runningKeyRef.current = pushRunning('run workflow');
    try {
      if (withSave && save) {
        await save();
      }
      const result = await startWorkflow({ projectName, workflowName });
      if (result.status === 'busy') {
        const message = 'The server is busy running another workflow. Try again shortly.';
        setError(message);
        pushError(`run workflow: ${message}`);
        finishRun();
      } else if (result.token) {
        // The poll effect takes over from here.
        setToken(result.token);
      } else {
        finishRun();
      }
    } catch (e: any) {
      setError(e?.message ?? String(e));
      pushError(`run workflow: ${e?.message ?? e}`);
      finishRun();
    }
  };

  // Only save on click when there is something to save.
  const saveOnClick = saveBeforeRun && canSave && Boolean(isChanged);

  // With unsaved changes and saving off, a plain run would execute the stale
  // saved version, so block it. The right-click menu still opens (it lives on
  // the wrapper Box) so the user can save-and-run or turn saving back on.
  const runBlocked = canSave && Boolean(isChanged) && !saveBeforeRun;

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

  const effectiveDisabled = disabled || runBlocked;
  let title = 'Run workflow (right click for more options)';
  if (disabled && disabledReason) {
    title = disabledReason;
  } else if (runBlocked) {
    title = 'Save changes before running (right click for options)';
  }
  const icon = saveBeforeRun && canSave
    ? (
      <Box sx={{ position: 'relative', display: 'inline-flex' }}>
        <PlayArrow />
        <Save sx={{ position: 'absolute', right: -5, bottom: -5, fontSize: 12 }} />
      </Box>
    )
    : <PlayArrow />;

  return (
    <>
      {/* The context menu lives on the wrapper so right click still opens it
          when the button itself is disabled by unsaved changes. */}
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
        <MenuItem onClick={() => runFromMenu(false)} disabled={canSave && Boolean(isChanged)}>
          Run
        </MenuItem>
        {canSave && (
          <MenuItem onClick={() => runFromMenu(true)} disabled={!isChanged}>
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
      <WorkflowOutputDialog
        open={open}
        running={running}
        output={output}
        error={error}
        workflowName={workflowName}
        onClose={() => { return setOpen(false); }}
      />
    </>
  );
};
