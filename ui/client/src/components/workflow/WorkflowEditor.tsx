import { Add } from '@mui/icons-material';
import { Box, Button, Stack, TextField, Typography } from '@mui/material';
import { useState } from 'react';
import { WorkflowBlock, WorkflowData, WorkflowNode } from '../../shared/types';
import { getWorkflowBlock, isTopLevelBlock } from '../../shared/workflow';
import { WorkflowGraph } from './WorkflowGraph';
import { WorkflowNodeEditor } from './WorkflowNodeEditor';

// Editor for a Hermes workflow document. Shows the solver and a node graph
// (edges from each node's `requires`); selecting a node opens it for editing.
export const WorkflowEditor = ({
  workflow,
  setWorkflow,
}: {
  workflow?: WorkflowData,
  setWorkflow: (workflow: WorkflowData) => void,
}) => {
  const [newNodeName, setNewNodeName] = useState('');
  const [selectedNode, setSelectedNode] = useState<string | undefined>(undefined);
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
    if (selectedNode === name) {
      setSelectedNode(undefined);
    }
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
            <WorkflowGraph
              nodeNames={nodeNames}
              nodes={block.nodes ?? {}}
              selectedNode={selectedNode}
              onSelectNode={setSelectedNode}
            />
            <Stack direction="row" spacing={1} alignItems="center" sx={{ mb: 2 }}>
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
            {selectedNode && block.nodes?.[selectedNode] && (
              <WorkflowNodeEditor
                name={selectedNode}
                node={block.nodes[selectedNode]}
                otherNodeNames={nodeNames.filter(n => n !== selectedNode)}
                setNode={(node) => setNode(selectedNode, node)}
                deleteNode={() => deleteNode(selectedNode)}
              />
            )}
          </>
        )
      }
    </Box>
  );
};
