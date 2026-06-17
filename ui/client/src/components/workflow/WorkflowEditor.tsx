import { Add } from '@mui/icons-material';
import { Box, Button, Stack, TextField, Typography } from '@mui/material';
import { useState } from 'react';
import { WorkflowBlock, WorkflowData, WorkflowNode } from '../../shared/types';
import { getWorkflowBlock, isTopLevelBlock } from '../../shared/workflow';
import { WorkflowNodeEditor } from './WorkflowNodeEditor';

// v1 editor for a Hermes workflow document. It shows the solver and lists the
// nodes in execution order; each node exposes its type and an editable tree of
// input parameters. The node-graph (ReactFlow) view will be a later iteration
// over this same data.
export const WorkflowEditor = ({
  workflow,
  setWorkflow,
}: {
  workflow?: WorkflowData,
  setWorkflow: (workflow: WorkflowData) => void,
}) => {
  const [newNodeName, setNewNodeName] = useState('');
  const block = getWorkflowBlock(workflow);

  // Preserve the original nesting (the block may be wrapped in an extra
  // { workflow: {...} } level) when writing edits back.
  const setBlock = (newBlock: WorkflowBlock) => {
    setWorkflow(isTopLevelBlock(workflow) ? newBlock : { ...workflow, workflow: newBlock });
  };

  const setNode = (name: string, node: WorkflowNode) => {
    setBlock({ ...block, nodes: { ...block?.nodes, [name]: node } });
  };

  const nodeNames = block?.nodeList?.length
    ? block.nodeList
    : Object.keys(block?.nodes ?? {});

  const addNode = () => {
    const name = newNodeName.trim();
    if (!name || nodeNames.includes(name)) return;
    setBlock({
      ...block,
      nodeList: [...nodeNames, name],
      nodes: { ...block?.nodes, [name]: { type: '', Execution: { input_parameters: {} } } },
    });
    setNewNodeName('');
  };

  const deleteNode = (name: string) => {
    const nodes = { ...block?.nodes };
    delete nodes[name];
    setBlock({ ...block, nodeList: nodeNames.filter(n => n !== name), nodes });
  };

  return (
    <Box sx={{ maxWidth: 900 }}>
      {!block
        ? <Typography color="text.secondary">No workflow found in this document.</Typography>
        : (
          <>
            <TextField
              label="Solver"
              size="small"
              value={block.solver ?? ''}
              onChange={(e) => setBlock({ ...block, solver: e.target.value })}
              sx={{ mb: 2 }}
            />
            <Typography variant="h6" sx={{ mb: 1 }}>
              Nodes ({nodeNames.length})
            </Typography>
            {nodeNames.map((name) => (
              <WorkflowNodeEditor
                key={name}
                name={name}
                node={block.nodes?.[name] ?? {}}
                otherNodeNames={nodeNames.filter(n => n !== name)}
                setNode={(node) => setNode(name, node)}
                deleteNode={() => deleteNode(name)}
              />
            ))}
            <Stack direction="row" spacing={1} alignItems="center" sx={{ mt: 1 }}>
              <TextField
                label="New node name"
                size="small"
                value={newNodeName}
                onChange={(e) => setNewNodeName(e.target.value)}
                onKeyDown={(e) => { if (e.key === 'Enter') addNode(); }}
              />
              <Button
                startIcon={<Add />}
                onClick={addNode}
                disabled={!newNodeName.trim() || nodeNames.includes(newNodeName.trim())}
              >
                Add node
              </Button>
            </Stack>
          </>
        )
      }
    </Box>
  );
};
