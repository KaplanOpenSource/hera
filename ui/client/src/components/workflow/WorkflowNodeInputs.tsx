import { SimpleTreeView } from '@mui/x-tree-view';
import { DetailsViewItem } from '../details/DetailsViewItem';
import { FieldDef } from '../details/fieldDef';

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
        // The chevron centers on a row that reserves extra space below for the
        // "required" helper, so it sits slightly low; nudge it up to the title.
        '& .MuiTreeItem-iconContainer': { transform: 'translateY(-4px)' },
      }}
    >
      <DetailsViewItem
        itemKey="input_parameters"
        itemValue={params}
        parentKey={undefined}
        def={paramsDef}
        setItemValue={onChangeParams}
      />
    </SimpleTreeView>
  );
};
