import { create } from 'zustand';

// A request to focus one node of one workflow. The counter makes a repeat pick of
// the same node a fresh request, so the canvas centres on it again.
export type NodeFocus = {
  workflowName: string,
  nodeName: string,
  seq: number,
};

type WorkflowFocusStore = {
  focus: NodeFocus | null,
  focusNode: (workflowName: string, nodeName: string) => void,
};

// Lets the output tab point the canvas tab at a node; the two are separate tabs,
// so they share this instead of props.
export const useWorkflowFocusStore = create<WorkflowFocusStore>((set) => {
  return {
    focus: null,
    focusNode: (workflowName, nodeName) => {
      return set((state) => {
        return { focus: { workflowName, nodeName, seq: (state.focus?.seq ?? 0) + 1 } };
      });
    },
  };
});
