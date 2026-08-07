import { WorkflowBlock, WorkflowData } from '../../shared/types';
import { getWorkflowBlock } from '../../shared/workflow';

// A named starter workflow the user can drop onto the canvas from the templates
// menu. Picking one replaces the current workflow's nodes with the block below.
export interface WorkflowTemplate {
  id: string;
  label: string;
  block: WorkflowBlock;
}

// Each JSON file in ./templates is a Hermes workflow. Drop a new file in that
// folder and it shows up in the templates menu — no code change needed. The
// file name (without .json) becomes the menu label. Keep them safe: nothing
// that runs a real simulation or touches the user's machine.
const templateFiles = import.meta.glob<WorkflowData>('./templates/*.json', {
  eager: true,
  import: 'default',
});

const labelFromPath = (path: string): string => {
  const file = path.split('/').pop() ?? path;
  return file.replace(/\.json$/, '');
};

export const workflowTemplates: WorkflowTemplate[] = Object.keys(templateFiles)
  .sort()
  .map(path => {
    return { label: labelFromPath(path), block: getWorkflowBlock(templateFiles[path]) };
  })
  .filter((t): t is { label: string, block: WorkflowBlock } => t.block !== undefined)
  .map(t => {
    return { id: t.label, label: t.label, block: t.block };
  });
