import { PlayArrow, Save } from '@mui/icons-material';
import { Box, Menu, MenuItem, SxProps, Theme } from '@mui/material';
import { MouseEvent, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { runWorkflow } from '../../io/runWorkflow';
import { dismiss, pushError, pushInfo, pushRunning } from '../../io/snackbar';
import { useViewSettingsStore } from '../../stores/useViewSettingsStore';
import { WorkflowOutputDialog } from './log/WorkflowOutputDialog';

// Runs a saved workflow via the server's /run_workflow endpoint. The run is
// synchronous — the output dialog opens right away and shows a spinner until the
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
  const [menuAnchor, setMenuAnchor] = useState<{ x: number, y: number } | null>(null);
  const saveBeforeRun = useViewSettingsStore(state => state.viewSettings.alwaysSaveBeforeRun);
  const setViewSettings = useViewSettingsStore(state => state.setViewSettings);

  const canSave = Boolean(save);

  const doRun = async (withSave: boolean) => {
    setOpen(true);
    setRunning(true);
    setOutput(null);
    setError(null);
    const key = pushRunning('run workflow');
    try {
      if (withSave && save) {
        await save();
      }
      const { output } = await runWorkflow({ projectName, workflowName });
      setOutput(output ?? '');
      pushInfo(`Workflow "${workflowName}" finished`);
    } catch (e: any) {
      setError(e?.message ?? String(e));
      pushError(`run workflow: ${e?.message ?? e}`);
    } finally {
      setRunning(false);
      dismiss(key);
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
