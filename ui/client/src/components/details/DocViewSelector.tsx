import { AccountTree, DataObject, DynamicForm, Science } from '@mui/icons-material';
import { ToggleButton, ToggleButtonGroup, Tooltip } from '@mui/material';
import { TabKind } from '../../shared/tabKind';

// The view modes: the always-present presentation modes (raw, formulated) plus
// the document-kind editors, which reuse the TabKind enum.
export type DocView = 'raw' | 'formulated' | TabKind.Agent | TabKind.Workflow;

// MUI removes pointer events on disabled buttons, which also suppresses their
// tooltip. Re-enable them so the explanatory tooltip still shows on hover.
const DISABLED_TOOLTIP_SX = { '&.Mui-disabled': { pointerEvents: 'auto' } };

export const DocViewSelector = ({
  docView,
  setDocView,
  enabled,
}: {
  docView: DocView,
  setDocView: (docView: DocView) => void,
  enabled: Partial<Record<TabKind, boolean>>,
}) => {
  return (
    <ToggleButtonGroup
      exclusive
      color="primary"
      size="small"
      value={docView}
      onChange={(_e, v: DocView | null) => { if (v) setDocView(v); }}
      sx={{ mr: 1 }}
    >
      <ToggleButton value="raw">
        <Tooltip title="Raw — the document exactly as stored, with nothing hidden">
          <DataObject fontSize="small" />
        </Tooltip>
      </ToggleButton>
      <ToggleButton value="formulated">
        <Tooltip title="Formulated — readable view with boilerplate metadata hidden">
          <DynamicForm fontSize="small" />
        </Tooltip>
      </ToggleButton>
      <ToggleButton value={TabKind.Agent} disabled={!enabled[TabKind.Agent]} sx={DISABLED_TOOLTIP_SX}>
        <Tooltip title={enabled[TabKind.Agent]
          ? 'Agent — edit effects and physical properties'
          : 'Agent — only available for agent documents'}>
          <Science fontSize="small" />
        </Tooltip>
      </ToggleButton>
      <ToggleButton value={TabKind.Workflow} disabled={!enabled[TabKind.Workflow]} sx={DISABLED_TOOLTIP_SX}>
        <Tooltip title={enabled[TabKind.Workflow]
          ? 'Workflow — edit solver, nodes and parameters'
          : 'Workflow — only available for workflow documents'}>
          <AccountTree fontSize="small" />
        </Tooltip>
      </ToggleButton>
    </ToggleButtonGroup>
  );
};
