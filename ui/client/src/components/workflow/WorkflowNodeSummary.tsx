import { Divider, Stack, Typography, useTheme } from '@mui/material';
import { Handle, Position } from '@xyflow/react';
import { FieldDef } from '../details/fieldDef';
import { FieldSourceDot } from '../details/FieldSourceDot';
import { friendlyParamName } from './friendlyParamName';
import { INPUT } from './references/knownKinds';
import { inputHandleId, outputHandleId } from './workflowDataflow';
import { WorkflowRequiresLine } from './WorkflowRequiresLine';

// Puts a small dot on the card's right edge, where a dataflow line leaves from.
const SOURCE_DOT_STYLE = { position: 'relative', top: 'auto', right: -12, transform: 'none', width: 8, height: 8, border: 'none' } as const;
// Wraps the row's own source dot, so a dataflow line lands on it.
const TARGET_DOT_STYLE = { position: 'relative', top: 'auto', left: 'auto', transform: 'none', width: 'auto', height: 'auto', minWidth: 0, minHeight: 0, background: 'transparent', border: 'none', borderRadius: 0 } as const;

// The small card a node shows while the pointer is away from it: its name, its
// type, one row per parameter - the same dot and name as the editor, with no
// value - and, under a divider, one row per output. Hovering the node lays the
// full editor over this, so nothing here is editable.
export const WorkflowNodeSummary = ({
  name,
  type,
  params,
  paramsDef,
  requires,
  outputs,
  withHandles,
}: {
  name: string,
  type?: string,
  params: { [key: string]: unknown },
  paramsDef: FieldDef,
  requires: string[],
  outputs: string[],
  // False while the editor is open over the card, which then has the dots.
  withHandles: boolean,
}) => {
  const theme = useTheme();
  // The same line dots as the editor, so lines reach the rows while the node is small.
  const sourceDotStyle = { ...SOURCE_DOT_STYLE, background: theme.palette.primary.main };
  return (
    <Stack spacing={0.25} sx={{ minWidth: 0 }}>
      <Typography noWrap sx={{ fontSize: 13, fontWeight: 600 }}>
        {name}
      </Typography>
      <Typography noWrap variant="caption" color={type ? 'text.secondary' : 'warning.main'}>
        {type || 'no type'}
      </Typography>
      <WorkflowRequiresLine requires={requires} />
      {Object.keys(params).map(param => {
        const dot = <FieldSourceDot source={paramsDef.children?.[param]?.source} showUnknown />;
        // A dict or list holds no value of its own, so no new line may land on it.
        const value = params[param];
        const isLeaf = typeof value !== 'object' || value === null;
        return (
          <Stack key={param} direction="row" spacing={1} sx={{ alignItems: 'center', minWidth: 0 }}>
            {withHandles && (
              <Handle type="target" id={inputHandleId(name, param)} position={Position.Left} isConnectableEnd={isLeaf} style={TARGET_DOT_STYLE}>
                {dot}
              </Handle>
            )}
            {!withHandles && dot}
            <Typography noWrap variant="caption" color="text.secondary" sx={{ flexGrow: 1 }}>
              {friendlyParamName(param)}
            </Typography>
            {withHandles && (
              <Handle type="source" id={INPUT.handleId(name, param)} position={Position.Right} style={sourceDotStyle} />
            )}
          </Stack>
        );
      })}
      {outputs.length > 0 && <Divider sx={{ my: 0.25 }} />}
      {outputs.map(output => (
        <Stack key={output} direction="row" spacing={1} sx={{ alignItems: 'center', justifyContent: 'flex-end', minWidth: 0 }}>
          <Typography noWrap variant="caption" color="text.secondary">
            {output}
          </Typography>
          {withHandles && (
            <Handle type="source" id={outputHandleId(name, output)} position={Position.Right} style={sourceDotStyle} />
          )}
        </Stack>
      ))}
    </Stack>
  );
};
