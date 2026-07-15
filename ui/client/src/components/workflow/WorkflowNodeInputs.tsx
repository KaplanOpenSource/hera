import { SimpleTreeView } from '@mui/x-tree-view';
import { Handle, Position } from '@xyflow/react';
import { DetailsViewItem, keyForDetailsViewItem } from '../details/DetailsViewItem';
import { FieldDef } from '../details/fieldDef';
import { FieldSourceDot } from '../details/FieldSourceDot';
import { inputHandleId } from './workflowDataflow';

// The node field holding a workflow node's input parameters — the key of the tree
// this component renders, and the parent key of each top-level parameter row.
export const INPUT_PARAMETERS_KEY = 'input_parameters';

// The node's input_parameters, shown as an editable tree. Expansion is
// controlled by the parent so the same chevron can also show/hide the outputs.
export const WorkflowNodeInputs = ({
  params,
  paramsDef,
  expandedItems,
  onExpandedItemsChange,
  onChangeParams,
  onFieldContextMenu,
  onFieldInlineEdit,
}: {
  params: { [key: string]: any },
  paramsDef: FieldDef,
  expandedItems: string[],
  onExpandedItemsChange: (itemIds: string[]) => void,
  onChangeParams: (newParams: any) => void,
  // Right-click on a top-level parameter row: open a menu for that field. The
  // caret is the click's position within the input value, when it lands on one.
  onFieldContextMenu: (param: string, x: number, y: number, caret?: number) => void,
  // Typing / caret moves in a top-level parameter's value editor, for inline
  // reference autocomplete: the param, the current value, the caret position, and
  // the input element (to anchor the suggestion menu to).
  onFieldInlineEdit: (param: string, value: string, caret: number | null, el: HTMLInputElement) => void,
}) => {
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
        setItemValue={onChangeParams}
        // Right-click on a top-level parameter opens a menu for that field.
        // Stop the event so ReactFlow's node menu doesn't also open; other rows
        // (the title, nested keys) fall through to the node menu.
        onRowContextMenu={(itemKey, parentKey, event) => {
          if (parentKey === keyForDetailsViewItem(INPUT_PARAMETERS_KEY)) {
            event.preventDefault();
            event.stopPropagation();
            // When the right-click lands on the input itself, capture the caret
            // so a reference can be inserted at that spot in the value.
            const el = event.target as HTMLInputElement;
            const caret = typeof el.selectionStart === 'number' ? el.selectionStart : undefined;
            onFieldContextMenu(itemKey, event.clientX, event.clientY, caret);
          }
        }}
        // Typing in a top-level parameter's value reports the caret so the editor
        // can offer inline reference suggestions; nested rows are ignored.
        onValueCaret={(itemKey, parentKey, value, caret, el) => {
          if (parentKey === keyForDetailsViewItem(INPUT_PARAMETERS_KEY)) {
            onFieldInlineEdit(itemKey, value, caret, el);
          }
        }}
        // Each parameter row shows its source dot before the name; on top-level
        // rows wrap that dot in a target handle (id = the parameter name) so a
        // dataflow line from another node's output can land on it.
        renderBeforeName={(itemKey, parentKey, def) => (
          parentKey === keyForDetailsViewItem(INPUT_PARAMETERS_KEY) ? (
            <Handle
              type="target"
              id={inputHandleId(itemKey)}
              position={Position.Left}
              style={{ position: 'relative', top: 'auto', left: 'auto', transform: 'none', width: 'auto', height: 'auto', minWidth: 0, minHeight: 0, background: 'transparent', border: 'none', borderRadius: 0 }}
            >
              <FieldSourceDot source={def?.source} showUnknown />
            </Handle>
          ) : (
            <FieldSourceDot source={def?.source} />
          )
        )}
      />
    </SimpleTreeView>
  );
};
