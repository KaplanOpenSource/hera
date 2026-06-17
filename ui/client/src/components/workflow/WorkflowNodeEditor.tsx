import { Delete } from '@mui/icons-material';
import { Box, Divider, Stack, TextField, Typography } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { WorkflowNode } from '../../shared/types';
import { DetailsViewItem, keyForDetailsViewItem } from '../details/DetailsViewItem';

// Editor for a single workflow node: its type and an editable tree of its input
// parameters, plus a delete control. (Requires edges are drawn on the canvas.)
export const WorkflowNodeEditor = ({
  name,
  node,
  setNode,
  deleteNode,
}: {
  name: string,
  node: WorkflowNode,
  setNode: (node: WorkflowNode) => void,
  deleteNode: () => void,
}) => {
  const params = node.Execution?.input_parameters ?? {};

  return (
    <Box sx={{ mb: 2 }}>
      <Stack direction="row" spacing={2} alignItems="center">
        <Typography variant="subtitle1" sx={{ fontWeight: 600 }}>
          {name}
        </Typography>
        <TextField
          label="type"
          size="small"
          value={node.type ?? ''}
          onChange={(e) => setNode({ ...node, type: e.target.value })}
        />
        <ButtonTooltip
          title={'Delete ' + name}
          onClick={deleteNode}
        >
          <Delete />
        </ButtonTooltip>
      </Stack>
      <SimpleTreeView
        defaultExpandedItems={[keyForDetailsViewItem('input_parameters')]}
      >
        <DetailsViewItem
          itemKey="input_parameters"
          itemValue={params}
          parentKey={undefined}
          setItemValue={(newVal) => setNode({
            ...node,
            Execution: { ...node.Execution, input_parameters: newVal },
          })}
        />
      </SimpleTreeView>
      <Divider sx={{ mt: 1 }} />
    </Box>
  );
};
