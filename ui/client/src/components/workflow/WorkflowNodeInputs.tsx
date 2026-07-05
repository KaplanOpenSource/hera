import { SimpleTreeView } from '@mui/x-tree-view';
import { Handle, Position } from '@xyflow/react';
import { DetailsViewItem, keyForDetailsViewItem } from '../details/DetailsViewItem';
import { FieldDef } from '../details/fieldDef';
import { FieldSourceDot } from '../details/FieldSourceDot';
import { inputHandleId } from './workflowDataflow';

// The node's input_parameters, shown as an editable tree. Expansion is
// controlled by the parent so the same chevron can also show/hide the outputs.
export const WorkflowNodeInputs = ({
  params,
  paramsDef,
  expandedItems,
  onExpandedItemsChange,
  onChangeParams,
}: {
  params: { [key: string]: any },
  paramsDef: FieldDef,
  expandedItems: string[],
  onExpandedItemsChange: (itemIds: string[]) => void,
  onChangeParams: (newParams: any) => void,
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
        itemKey="input_parameters"
        itemValue={params}
        parentKey={undefined}
        def={paramsDef}
        setItemValue={onChangeParams}
        // Each parameter row shows its source dot before the name; on top-level
        // rows wrap that dot in a target handle (id = the parameter name) so a
        // dataflow line from another node's output can land on it.
        renderBeforeName={(itemKey, parentKey, def) => (
          parentKey === keyForDetailsViewItem('input_parameters') ? (
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
