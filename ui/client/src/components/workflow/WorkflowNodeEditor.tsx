import { Delete } from '@mui/icons-material';
import { Box, Divider, Stack, TextField, Typography } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { WorkflowNode } from '../../shared/types';
import { DetailsViewItem, keyForDetailsViewItem } from '../details/DetailsViewItem';
import { WorkflowRequiresSelect } from './WorkflowRequiresSelect';

// Editor for a single workflow node: its type, the nodes it requires, and an
// editable tree of its input parameters, plus a delete control.
export const WorkflowNodeEditor = ({
  name,
  node,
  otherNodeNames,
  setNode,
  deleteNode,
}: {
  name: string,
  node: WorkflowNode,
  otherNodeNames: string[],
  setNode: (node: WorkflowNode) => void,
  deleteNode: () => void,
}) => {
  const params = node.Execution?.input_parameters ?? {};

  // Drop the key entirely when nothing is required, otherwise store the list.
  const setRequires = (values: string[]) => {
    const { requires: _drop, ...rest } = node;
    setNode(values.length ? { ...rest, requires: values } : rest);
  };

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
        <WorkflowRequiresSelect
          requires={node.requires}
          options={otherNodeNames}
          onChange={setRequires}
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
