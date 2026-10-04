import { Typography, useTheme } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { Handle, Position } from '@xyflow/react';
import { DetailsViewItem, keyForDetailsViewItem } from '../details/DetailsViewItem';
import { FieldDef } from '../details/fieldDef';
import { FieldSourceDot } from '../details/FieldSourceDot';
import { INPUT } from './references/knownKinds';
import { inputHandleId } from './workflowDataflow';
import { paramPathOf, valueAtPath } from './paramPath';
import { friendlyParamName } from './friendlyParamName';

// The node field holding a workflow node's input parameters — the key of the tree
// this component renders, and the parent key of each top-level parameter row.
export const INPUT_PARAMETERS_KEY = 'input_parameters';

// The tree key of the section every parameter row sits under.
const INPUTS_KEY = keyForDetailsViewItem(INPUT_PARAMETERS_KEY);

// Where a row's value sits inside input_parameters, or null when the row is not
// one of its values (the section's own title row). A row nested in a dict or a
// list gets the full dotted path, so a line can land on it too.
export const paramPathOfRow = (itemKey: string, parentKey: string | undefined): string | null => {
  if (parentKey === INPUTS_KEY) {
    return itemKey;
  }
  if (parentKey?.startsWith(`${INPUTS_KEY}/`)) {
    return paramPathOf([...parentKey.slice(INPUTS_KEY.length + 1).split('/'), itemKey]);
  }
  return null;
};

// A row that holds a value of its own, rather than a dict or a list of more
// rows. Only a leaf can be the target of a dataflow line.
const isLeafValue = (value: any): boolean => {
  return value === null || typeof value !== 'object';
};

// The node's input_parameters, shown as an editable tree. Expansion is
// controlled by the parent so the same chevron can also show/hide the outputs.
export const WorkflowNodeInputs = ({
  nodeName,
  params,
  paramsDef,
  expandedItems,
  onExpandedItemsChange,
  onChangeParams,
  onFieldContextMenu,
  onFieldInlineEdit,
}: {
  nodeName: string,
  params: { [key: string]: any },
  paramsDef: FieldDef,
  expandedItems: string[],
  onExpandedItemsChange: (itemIds: string[]) => void,
  onChangeParams: (newParams: any) => void,
  // Right-click on a parameter row: open a menu for that field, named by its
  // path. The caret is the click's position within the input value, when it
  // lands on one.
  onFieldContextMenu: (paramPath: string, x: number, y: number, caret?: number) => void,
  // Typing / caret moves in a parameter's value editor, for inline reference
  // autocomplete: the row's path, the current value, the caret position, and
  // the input element (to anchor the suggestion menu to).
  onFieldInlineEdit: (paramPath: string, value: string, caret: number | null, el: HTMLInputElement) => void,
}) => {
  const theme = useTheme();
  // The path of a row that can take a dataflow line, or null for a row that
  // can't: the section title, or a dict / list row that holds no value.
  const leafPathOfRow = (itemKey: string, parentKey: string | undefined): string | null => {
    const path = paramPathOfRow(itemKey, parentKey);
    if (path === null || !isLeafValue(valueAtPath(params, path))) {
      return null;
    }
    return path;
  };
  return (
    <SimpleTreeView
      expandedItems={expandedItems}
      onExpandedItemsChange={(_e, itemIds) => onExpandedItemsChange(itemIds)}
      sx={{
        flexGrow: 1,
        minWidth: 0,
        '& .MuiTreeItem-label .MuiTypography-root': { fontSize: '0.875rem' },
        // The title row (input_parameters) inherits DetailsViewItem's tall row
        // margins, which exist to reserve space under leaf fields; the title has
        // no field, so strip them on the top-level row to keep it compact.
        '& > .MuiTreeItem-root > .MuiTreeItem-content .MuiTreeItem-label > .MuiStack-root': {
          marginTop: '0 !important',
          marginBottom: '0 !important',
        },
        // The chevron centers on a row that reserves extra space below for the
        // "required" helper, so it sits slightly low; nudge it up to the title.
        '& .MuiTreeItem-iconContainer': { transform: 'translateY(-2px)' },
      }}
    >
      <DetailsViewItem
        itemKey={INPUT_PARAMETERS_KEY}
        itemValue={params}
        parentKey={undefined}
        def={paramsDef}
        // The title row only names the section, so no hermes name and no type chip.
        // A parameter row shows a readable name; its hermes name shows while editing.
        nameForView={(itemKey, parentKey) => {
          if (parentKey === undefined) {
            return <Typography sx={{ whiteSpace: 'nowrap', flexShrink: 0 }}>Parameters</Typography>;
          }
          if (parentKey !== INPUTS_KEY) {
            return undefined;
          }
          return (
            <Typography sx={{ whiteSpace: 'nowrap', minWidth: '100px', flexShrink: 0 }}>
              {friendlyParamName(itemKey)}
            </Typography>
          );
        }}
        hideTypeSelector
        setItemValue={onChangeParams}
        // Right-click on any parameter row opens a menu for that field. Stop the
        // event so ReactFlow's node menu doesn't also open; the section title
        // row falls through to the node menu.
        onRowContextMenu={(itemKey, parentKey, event) => {
          const path = paramPathOfRow(itemKey, parentKey);
          if (path !== null) {
            event.preventDefault();
            event.stopPropagation();
            // When the right-click lands on the input itself, capture the caret
            // so a reference can be inserted at that spot in the value.
            const el = event.target as HTMLInputElement;
            const caret = typeof el.selectionStart === 'number' ? el.selectionStart : undefined;
            onFieldContextMenu(path, event.clientX, event.clientY, caret);
          }
        }}
        // Typing in a parameter's value reports the caret so the editor can
        // offer inline reference suggestions.
        onValueCaret={(itemKey, parentKey, value, caret, el) => {
          const path = paramPathOfRow(itemKey, parentKey);
          if (path !== null) {
            onFieldInlineEdit(path, value, caret, el);
          }
        }}
        // Each parameter row shows its source dot before the name; on a leaf row
        // wrap that dot in a target handle (id = the row's path) so a dataflow
        // line from another node's output can land on it, however deeply the
        // value sits. A dict or list row holds no value of its own, so it gets
        // a plain dot.
        renderBeforeName={(itemKey, parentKey, def) => {
          const path = leafPathOfRow(itemKey, parentKey);
          if (path === null) {
            return <FieldSourceDot source={def?.source} />;
          }
          return (
            <Handle
              type="target"
              id={inputHandleId(nodeName, path)}
              position={Position.Left}
              style={{ position: 'relative', top: 'auto', left: 'auto', transform: 'none', width: 'auto', height: 'auto', minWidth: 0, minHeight: 0, background: 'transparent', border: 'none', borderRadius: 0 }}
            >
              <FieldSourceDot source={def?.source} showUnknown />
            </Handle>
          );
        }}
        // A top-level parameter also gets a source dot on the node's right edge,
        // so another node can read this parameter's value. Only top-level: a
        // reference can only name a parameter, not a value inside one.
        renderAfterValue={(itemKey, parentKey) => {
          if (parentKey !== INPUTS_KEY) {
            return undefined;
          }
          return (
            <Handle
              type="source"
              id={INPUT.handleId(nodeName, itemKey)}
              position={Position.Right}
              style={{ position: 'relative', top: 'auto', right: -14, transform: 'none', width: 8, height: 8, background: theme.palette.primary.main, border: 'none' }}
            />
          );
        }}
      />
    </SimpleTreeView>
  );
};
