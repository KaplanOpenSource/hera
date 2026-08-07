import { Box, Typography } from '@mui/material';
import { ReactNode, useState } from 'react';
import { WorkflowBlock, WorkflowData, WorkflowNode } from '../../shared/types';
import { getWorkflowBlock, isTopLevelBlock, normalizeRequires } from '../../shared/workflow';
import { NodeCatalogReader, useNodeCatalog } from './useNodeCatalog';
import { WorkflowGraph } from './WorkflowGraph';

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
  actionButtons,
}: {
  workflow?: WorkflowData,
  setWorkflow: (workflow: WorkflowData) => void,
  actionButtons?: ReactNode,
}) => {
  const [selectedNode, setSelectedNode] = useState<string | undefined>(undefined);
  const catalog = useNodeCatalog(s => s.catalog);
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

  // Replace the whole workflow with a starter template's block.
  const applyTemplate = (templateBlock: WorkflowBlock) => {
    setBlock(templateBlock);
    setSelectedNode(undefined);
  };

  const deleteNode = (name: string) => {
    const nodes = { ...block?.nodes };
    delete nodes[name];
    setBlock({ ...block, nodeList: nodeNames.filter(n => n !== name), nodes });
    if (selectedNode === name) {
      setSelectedNode(undefined);
    }
  };

  // Add/remove a requires edge (source must run before target).
  const addRequire = (source: string, target: string) => {
    const node = block?.nodes?.[target];
    if (!node) {
      return;
    }
    const current = normalizeRequires(node.requires);
    if (current.includes(source)) {
      return;
    }
    setNode(target, { ...node, requires: [...current, source] });
  };

  const removeRequire = (source: string, target: string) => {
    const node = block?.nodes?.[target];
    if (!node) {
      return;
    }
    const requires = normalizeRequires(node.requires).filter(r => r !== source);
    const { requires: _drop, ...rest } = node;
    setNode(target, requires.length > 0 ? { ...rest, requires } : rest);
  };

  return (
    <Box sx={{ flex: 1, display: 'flex', flexDirection: 'column', minHeight: 0 }}>
      <NodeCatalogReader />
      {!block
        ? <Typography color="text.secondary">No workflow found in this document.</Typography>
        : (
          <>
            <WorkflowGraph
              catalog={catalog}
              nodeNames={nodeNames}
              nodes={block.nodes ?? {}}
              selectedNode={selectedNode}
              actionButtons={actionButtons}
              onSelectNode={setSelectedNode}
              onAddNode={addNode}
              onApplyTemplate={applyTemplate}
              onRenameNode={renameNode}
              onSetNode={setNode}
              onAddRequire={addRequire}
              onRemoveRequire={removeRequire}
              onDeleteNode={deleteNode}
            />
          </>
        )
      }
    </Box>
  );
};
