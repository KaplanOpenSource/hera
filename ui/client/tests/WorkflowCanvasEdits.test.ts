import { describe, it, expect, vi } from 'vitest';
import { WorkflowCanvasEdits } from '../src/components/workflow/WorkflowCanvasEdits';
import { inputHandleId, nodeInputHandleId, nodeOutputHandleId, outputHandleId } from '../src/components/workflow/workflowDataflow';

const nodes = {
  A: { type: 'maker' },
  B: { type: 'plain', Execution: { input_parameters: { cmd: 'run {A.output.result} now', keep: 'as is' } } },
};

const build = () => {
  const callbacks = { onSetNode: vi.fn(), onAddRequire: vi.fn(), onRemoveRequire: vi.fn() };
  return { callbacks, edits: new WorkflowCanvasEdits(['A', 'B'], nodes, callbacks) };
};

// The parameters the last onSetNode call wrote.
const writtenParams = (onSetNode: ReturnType<typeof vi.fn>) => {
  return onSetNode.mock.calls.at(-1)![1].Execution.input_parameters;
};

describe('WorkflowCanvasEdits', () => {
  it('adds a requires link for a node-to-node line', () => {
    const { callbacks, edits } = build();
    edits.connect({ source: 'A', target: 'B', sourceHandle: nodeOutputHandleId('A'), targetHandle: nodeInputHandleId('B') });
    expect(callbacks.onAddRequire).toHaveBeenCalledWith('A', 'B');
  });

  it('writes a reference for an output-to-input line', () => {
    const { callbacks, edits } = build();
    edits.connect({ source: 'A', target: 'B', sourceHandle: outputHandleId('A', 'result'), targetHandle: inputHandleId('B', 'cmd') });
    expect(callbacks.onAddRequire).not.toHaveBeenCalled();
    expect(writtenParams(callbacks.onSetNode).cmd).toContain('{A.output.result}');
  });

  it('ignores a line with no ends', () => {
    const { callbacks, edits } = build();
    edits.connect({ source: null, target: null, sourceHandle: null, targetHandle: null });
    expect(callbacks.onAddRequire).not.toHaveBeenCalled();
    expect(callbacks.onSetNode).not.toHaveBeenCalled();
  });

  it('allows an output-to-input line without the cycle check', () => {
    const { edits } = build();
    expect(edits.canConnect({
      source: 'B',
      target: 'A',
      sourceHandle: outputHandleId('B', 'out'),
      targetHandle: inputHandleId('A', 'p'),
    })).toBe(true);
  });

  it('clears the reference a dataflow line stands for', () => {
    const { callbacks, edits } = build();
    edits.removeDataflowEdge('df:A.result->B.cmd');
    expect(writtenParams(callbacks.onSetNode).cmd).not.toContain('{A.output.result}');
  });

  it('removes a requires line by its ends', () => {
    const { callbacks, edits } = build();
    edits.removeEdges([{ id: 'A->B', source: 'A', target: 'B' }]);
    expect(callbacks.onRemoveRequire).toHaveBeenCalledWith('A', 'B');
  });

  it('deletes one field of a node', () => {
    const { callbacks, edits } = build();
    edits.deleteField('B', 'cmd');
    expect(Object.keys(writtenParams(callbacks.onSetNode))).toEqual(['keep']);
  });

  it('inserts a reference at the caret', () => {
    const { callbacks, edits } = build();
    edits.referenceOutput('B', 'keep', 'A', 'result', 2);
    expect(writtenParams(callbacks.onSetNode).keep).toBe('as{A.output.result} is');
  });
});
