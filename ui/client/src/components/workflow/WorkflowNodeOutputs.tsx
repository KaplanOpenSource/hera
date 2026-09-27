import { Typography } from '@mui/material';
import { SimpleTreeView, TreeItem } from '@mui/x-tree-view';
import { useState } from 'react';
import { keyForDetailsViewItem } from '../details/DetailsViewItem';
import { WorkflowNodeOutputChip } from './WorkflowNodeOutputChip';

// The key of the outputs tree, and the parent key of each output row.
export const OUTPUTS_KEY = 'outputs';

// The node's outputs, stacked under the inputs in their own collapsible list.
// One row per output, each with a single source dot on the node's right edge.
export const WorkflowNodeOutputs = ({
  nodeName,
  outputs,
}: {
  nodeName: string,
  outputs: string[],
}) => {
  const [expandedItems, setExpandedItems] = useState<string[]>([OUTPUTS_KEY]);
  return (
    <SimpleTreeView
      className="nodrag"
      expandedItems={expandedItems}
      onExpandedItemsChange={(_e, itemIds) => setExpandedItems(itemIds)}
      sx={{
        flexGrow: 1,
        minWidth: 0,
        '& .MuiTreeItem-label .MuiTypography-root': { fontSize: '0.875rem' },
        // MUI clips the label, which would cut the dot hanging off the right.
        '& .MuiTreeItem-content > .MuiTreeItem-label': { overflow: 'visible' },
      }}
    >
      <TreeItem
        itemId={OUTPUTS_KEY}
        label={<Typography sx={{ whiteSpace: 'nowrap' }}>Outputs</Typography>}
      >
        {outputs.map(name => (
          <TreeItem
            key={name}
            itemId={keyForDetailsViewItem(name, OUTPUTS_KEY)}
            label={<WorkflowNodeOutputChip nodeName={nodeName} name={name} />}
          />
        ))}
      </TreeItem>
    </SimpleTreeView>
  );
};
