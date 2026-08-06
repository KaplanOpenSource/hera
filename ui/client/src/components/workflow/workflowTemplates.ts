import { WorkflowBlock } from '../../shared/types';

// A named starter workflow the user can drop onto the canvas from the templates
// menu. Picking one replaces the current workflow's nodes with the block below.
export interface WorkflowTemplate {
  id: string;
  label: string;
  block: WorkflowBlock;
}

// Safe templates only — nothing that runs a real simulation or touches the
// user's machine. Good for smoke-testing "run workflow from the UI".
export const workflowTemplates: WorkflowTemplate[] = [
  {
    id: 'helloWorld',
    label: 'Hello World',
    block: {
      nodeList: ['SayHello'],
      nodes: {
        SayHello: {
          type: 'general.RunOsCommand',
          Execution: {
            input_parameters: {
              Method: 'Commands list',
              Command: ["echo 'hello from hera'", 'date'],
            },
          },
        },
      },
    },
  },
];
