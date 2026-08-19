// Helpers for Hermes workflow documents.
//
// A Hermes workflow is stored inside a Simulations document under
// `desc.workflow`. That value sometimes wraps the actual block in an extra
// `{ workflow: {...} }` level. The block of interest holds `nodeList` (the
// ordered node names) and `nodes` (each node's type + input parameters).

import { WorkflowBlock, WorkflowData, WorkflowDesc } from './types';

// Document type string set by the hermes workflow toolkit (DOCTYPE_WORKFLOW).
export const WORKFLOW_DOC_TYPE = 'hermesWorkflow';

// True for a parameter key that names a project: `ProjectName`, `projectName`,
// `projectname` - any casing. Hera/Hermes nodes use `ProjectName`.
export const isProjectNameKey = (key: string): boolean => {
  return /^projectname$/i.test(key);
};

// Sets every project-name parameter (any casing) in the block to the given
// project, so a template runs against the project it was dropped into. Looks at
// each node's top-level input_parameters and overwrites any existing value.
export const fillProjectName = (block: WorkflowBlock, projectName: string): WorkflowBlock => {
  const nodes = { ...block.nodes };
  for (const [name, node] of Object.entries(nodes)) {
    const params = node.Execution?.input_parameters;
    if (!params) {
      continue;
    }
    const nextParams = { ...params };
    let changed = false;
    for (const key of Object.keys(nextParams)) {
      if (isProjectNameKey(key)) {
        nextParams[key] = projectName;
        changed = true;
      }
    }
    if (changed) {
      nodes[name] = {
        ...node,
        Execution: { ...node.Execution, input_parameters: nextParams },
      };
    }
  }
  return { ...block, nodes };
};

// After a node's input_parameters were edited, set any project-name parameter
// that was just created or renamed into existence (a matching key not present
// before the edit) to the given project. A key that already existed is left
// alone, so the user can still override its value by hand.
export const fillNewProjectNameParams = (
  nextParams: { [key: string]: any },
  prevParams: { [key: string]: any },
  projectName: string,
): { [key: string]: any } => {
  let result = nextParams;
  for (const key of Object.keys(nextParams)) {
    if (isProjectNameKey(key) && !(key in prevParams)) {
      if (result === nextParams) {
        result = { ...nextParams };
      }
      result[key] = projectName;
    }
  }
  return result;
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
