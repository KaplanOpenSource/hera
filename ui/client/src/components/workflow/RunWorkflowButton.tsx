import { PlayArrow } from '@mui/icons-material';
import { SxProps, Theme } from '@mui/material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { fetchPython } from '../../io/fetchPython';
import { pushInfo } from '../../io/snackbar';
import { buildRunWorkflowCode } from './buildRunWorkflowCode';

// Runs a saved workflow by building and executing it from the DB. The run is
// synchronous — the button stays busy until the workflow finishes (errors are
// surfaced by fetchPython's snackbar).
export const RunWorkflowButton = ({
  projectName,
  workflowName,
  disabled,
  disabledReason,
  sx,
}: {
  projectName: string,
  workflowName: string,
  disabled?: boolean,
  disabledReason?: string,
  sx?: SxProps<Theme>,
}) => {
  const runWorkflow = async () => {
    const { data } = await fetchPython({
      results: ['dispatch_id'],
      label: 'run workflow',
      code: buildRunWorkflowCode({ projectName, workflowName }),
    });
    if (!data) {
      return;
    }
    pushInfo(`Workflow "${workflowName}" finished`);
  };

  const title = disabled && disabledReason ? disabledReason : 'Run workflow';
  return (
    <ButtonTooltip
      title={title}
      aria-label={title}
      disabled={disabled}
      onClick={runWorkflow}
      sx={sx}
    >
      <PlayArrow />
    </ButtonTooltip>
  );
};
