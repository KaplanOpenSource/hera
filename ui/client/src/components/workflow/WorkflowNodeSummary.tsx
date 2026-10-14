import { Stack, Typography } from '@mui/material';
import { FieldDef } from '../details/fieldDef';
import { FieldSourceDot } from '../details/FieldSourceDot';
import { friendlyParamName } from './friendlyParamName';
import { WorkflowRequiresLine } from './WorkflowRequiresLine';

// The small card a node shows while the pointer is away from it: its name, its
// type, and one row per parameter - the same dot and name as the editor, with
// no value. Hovering the node lays the full editor over this, so nothing here
// is editable.
export const WorkflowNodeSummary = ({
  name,
  type,
  paramNames,
  paramsDef,
  requires,
}: {
  name: string,
  type?: string,
  paramNames: string[],
  paramsDef: FieldDef,
  requires: string[],
}) => {
  return (
    <Stack spacing={0.25} sx={{ minWidth: 0 }}>
      <Typography noWrap sx={{ fontSize: 13, fontWeight: 600 }}>
        {name}
      </Typography>
      <Typography noWrap variant="caption" color={type ? 'text.secondary' : 'warning.main'}>
        {type || 'no type'}
      </Typography>
      <WorkflowRequiresLine requires={requires} />
      {paramNames.map(param => (
        <Stack key={param} direction="row" spacing={1} sx={{ alignItems: 'center', minWidth: 0 }}>
          <FieldSourceDot source={paramsDef.children?.[param]?.source} showUnknown />
          <Typography noWrap variant="caption" color="text.secondary">
            {friendlyParamName(param)}
          </Typography>
        </Stack>
      ))}
    </Stack>
  );
};
