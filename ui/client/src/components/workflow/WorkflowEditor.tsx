import { Box, TextField, Typography } from '@mui/material';
import { useState } from 'react';
import { WorkflowBlock, WorkflowData, WorkflowNode } from '../../shared/types';
import { getWorkflowBlock, isTopLevelBlock } from '../../shared/workflow';
import { WorkflowGraph } from './WorkflowGraph';
import { WorkflowNodeEditor } from './WorkflowNodeEditor';

// Returns the node's `requires` with oldName replaced by newName, preserving
// its single-name / list shape (or undefined when the node had no requires).
const normalizeRenamed = (requires: string | string[] | undefined, oldName: string, newName: string): string | string[] | undefined => {
  if (requires === undefined) {
    return undefined;
  }
  if (Array.isArray(requires)) {
    return requires.map(r => (r === oldName ? newName : r));
  }
  return requires === oldName ? newName : requires;
};

// Editor for a Hermes workflow document. Shows the solver and a node graph
// (edges from each node's `requires`); selecting a node opens it for editing.
export const WorkflowEditor = ({
  workflow,
  setWorkflow,
}: {
  workflow?: WorkflowData,
  setWorkflow: (workflow: WorkflowData) => void,
}) => {
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

  const uniqueNodeName = (): string => {
    let i = nodeNames.length + 1;
    let name = `node${i}`;
    while (nodeNames.includes(name)) {
      i += 1;
      name = `node${i}`;
    }
    return name;
  };

  const addNode = () => {
    const name = uniqueNodeName();
    setBlock({
      ...block,
      nodeList: [...nodeNames, name],
      nodes: { ...block?.nodes, [name]: { type: '', Execution: { input_parameters: {} } } },
    });
    setSelectedNode(name);
  };

  // Rename a node: update its key, the node list, and any `requires` references.
  const renameNode = (oldName: string, newName: string) => {
    const name = newName.trim();
    if (!name || name === oldName || nodeNames.includes(name)) {
      return;
    }
    const oldNodes = block?.nodes ?? {};
    const nodes: { [name: string]: WorkflowNode } = {};
    for (const key of Object.keys(oldNodes)) {
      const node = oldNodes[key];
      const requires = normalizeRenamed(node.requires, oldName, name);
      nodes[key === oldName ? name : key] = requires === undefined ? node : { ...node, requires };
    }
    setBlock({ ...block, nodeList: nodeNames.map(n => (n === oldName ? name : n)), nodes });
    if (selectedNode === oldName) {
      setSelectedNode(name);
    }
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
              onAddNode={addNode}
              onRenameNode={renameNode}
            />
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
