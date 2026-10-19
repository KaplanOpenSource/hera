import type { WorkflowChunk } from '../api/WorkflowApi';

export type { StartWorkflowResult, WorkflowChunk, PollWorkflowResult } from '../api/WorkflowApi';

// Joins per-task chunks into one flat text, in run order (for a plain-text view).
export const chunksToText = (
  chunks?: WorkflowChunk[] | null,
): string => {
  return (chunks ?? []).map((chunk) => { return chunk.text; }).join('');
};
