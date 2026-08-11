// Helpers for Hermes workflow documents.
//
// A Hermes workflow is stored inside a Simulations document under
// `desc.workflow`. That value sometimes wraps the actual block in an extra
// `{ workflow: {...} }` level. The block of interest holds `nodeList` (the
// ordered node names) and `nodes` (each node's type + input parameters).

import { WorkflowBlock, WorkflowData, WorkflowDesc } from './types';

// Document type string set by the hermes workflow toolkit (DOCTYPE_WORKFLOW).
export const WORKFLOW_DOC_TYPE = 'hermesWorkflow';

// Fills any node that takes a `projectName` parameter (but hasn't got one) with
// the given project, so a template runs against the project it was dropped into.
export const fillProjectName = (block: WorkflowBlock, projectName: string): WorkflowBlock => {
  const nodes = { ...block.nodes };
  for (const [name, node] of Object.entries(nodes)) {
    const params = node.Execution?.input_parameters?.Parameters;
    if (params && !params.projectName) {
      nodes[name] = {
        ...node,
        Execution: {
          ...node.Execution,
          input_parameters: {
            ...node.Execution?.input_parameters,
            Parameters: { ...params, projectName },
          },
        },
      };
    }
  }
  return { ...block, nodes };
};

// A node's `requires` as an array (it may be stored as a single name or a list).
export const normalizeRequires = (requires?: string | string[]): string[] => {
  if (requires === undefined) {
    return [];
  }
  return Array.isArray(requires) ? requires : [requires];
};

// Returns the inner workflow block (the one carrying nodeList/nodes),
// unwrapping the optional extra { workflow: {...} } level.
export const getWorkflowBlock = (workflow?: WorkflowData): WorkflowBlock | undefined => {
  if (!workflow) return undefined;
  if ('nodeList' in workflow || 'nodes' in workflow) return workflow;
  return getWorkflowBlock(workflow.workflow);
};

// True when the workflow value is the bare block (carries nodeList/nodes
// directly), false when the block is nested under an extra { workflow } level.
export const isTopLevelBlock = (workflow?: WorkflowData): boolean => {
  return !!workflow && ('nodeList' in workflow || 'nodes' in workflow);
};

// The solver name from the workflow block (empty string when unset).
export const getWorkflowSolver = (workflow?: WorkflowData): string => {
  return getWorkflowBlock(workflow)?.solver ?? '';
};

// Returns the workflow with its block's solver set, preserving the optional
// extra { workflow } nesting level.
export const setWorkflowSolver = (workflow: WorkflowData | undefined, solver: string): WorkflowData => {
  const block = { ...getWorkflowBlock(workflow), solver };
  return isTopLevelBlock(workflow) ? block : { ...workflow, workflow: block };
};

export const isWorkflowDoc = (doc?: { type?: string; desc?: WorkflowDesc }) => {
  if (!doc) return false;
  if (doc.type === WORKFLOW_DOC_TYPE) return true;
  return getWorkflowBlock(doc.desc?.workflow) !== undefined;
};
