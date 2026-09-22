import { describe, it, expect } from 'vitest';
import {
  deOverlappedFlowNodes,
  FitKind,
  flowMeasuredKey,
  flowNodeCenter,
  flowStructureKey,
  pendingFitAfterChange,
  rebuiltFlowNodes,
  restackedFlowNodes,
} from '../src/components/workflow/workflowFlowNodesUtils';

const flowNode = (id: string, x: number, y: number) => {
  return { id, type: 'workflow', position: { x, y }, data: {} };
};

describe('rebuiltFlowNodes', () => {
  it('keeps the position and size of a node that survives', () => {
    const prev = [{ ...flowNode('A', 5, 7), width: 300, height: 120 }];
    const built = rebuiltFlowNodes(prev, ['A'], { A: { x: 0, y: 0 } });
    expect(built[0].position).toEqual({ x: 5, y: 7 });
    expect(built[0].width).toBe(300);
    expect(built[0].height).toBe(120);
  });

  it('puts a new node where the layout says', () => {
    const built = rebuiltFlowNodes([], ['A'], { A: { x: 40, y: 80 } });
    expect(built[0].position).toEqual({ x: 40, y: 80 });
  });

  it('drops a node that is no longer in the workflow', () => {
    const built = rebuiltFlowNodes([flowNode('A', 0, 0), flowNode('B', 0, 0)], ['A'], { A: { x: 0, y: 0 } });
    expect(built.map((node) => { return node.id; })).toEqual(['A']);
  });
});

describe('restackedFlowNodes', () => {
  it('moves the nodes the layout knows and leaves the rest', () => {
    const moved = restackedFlowNodes([flowNode('A', 0, 0), flowNode('B', 1, 1)], { A: { x: 9, y: 9 } });
    expect(moved[0].position).toEqual({ x: 9, y: 9 });
    expect(moved[1].position).toEqual({ x: 1, y: 1 });
  });
});

describe('deOverlappedFlowNodes', () => {
  it('pushes a node down to its fixed y', () => {
    const fixed = deOverlappedFlowNodes([flowNode('A', 0, 0)], { A: { x: 0, y: 50 } });
    expect(fixed[0].position).toEqual({ x: 0, y: 50 });
  });

  it('returns the same list when nothing moved', () => {
    const prev = [flowNode('A', 0, 50)];
    expect(deOverlappedFlowNodes(prev, { A: { x: 0, y: 50 } })).toBe(prev);
  });
});

describe('pendingFitAfterChange', () => {
  it('fits everything on the first load', () => {
    expect(pendingFitAfterChange([], ['A'])).toEqual({ kind: FitKind.All });
  });

  it('pans to a single added node', () => {
    expect(pendingFitAfterChange(['A'], ['A', 'B'])).toEqual({ kind: FitKind.Node, nodeName: 'B' });
  });

  it('fits everything when several nodes arrive', () => {
    expect(pendingFitAfterChange(['A'], ['A', 'B', 'C'])).toEqual({ kind: FitKind.All });
  });

  it('fits everything when a node is swapped out', () => {
    expect(pendingFitAfterChange(['A'], ['B'])).toEqual({ kind: FitKind.All });
  });

  it('just refits when a node was removed', () => {
    expect(pendingFitAfterChange(['A', 'B'], ['A'])).toEqual({ kind: FitKind.Refit });
  });

  it('just refits when the names did not change', () => {
    expect(pendingFitAfterChange(['A'], ['A'])).toEqual({ kind: FitKind.Refit });
  });
});

describe('flowNodeCenter', () => {
  it('is the middle of the measured box', () => {
    expect(flowNodeCenter({ position: { x: 10, y: 20 }, measured: { width: 100, height: 40 } })).toEqual({ x: 60, y: 40 });
  });

  it('is the position itself before the node is measured', () => {
    expect(flowNodeCenter({ position: { x: 10, y: 20 } })).toEqual({ x: 10, y: 20 });
  });
});

describe('change keys', () => {
  it('changes when a node type changes', () => {
    const before = flowStructureKey(['A'], { A: { type: 'one' } }, []);
    const after = flowStructureKey(['A'], { A: { type: 'two' } }, []);
    expect(before).not.toBe(after);
  });

  it('does not change when only a parameter value changes', () => {
    const before = flowStructureKey(['A'], { A: { type: 'one', Execution: { input_parameters: { x: 1 } } } }, []);
    const after = flowStructureKey(['A'], { A: { type: 'one', Execution: { input_parameters: { x: 2 } } } }, []);
    expect(before).toBe(after);
  });

  it('rounds the measured heights', () => {
    const a = flowMeasuredKey([{ ...flowNode('A', 0, 0), measured: { height: 100.2 } }]);
    const b = flowMeasuredKey([{ ...flowNode('A', 0, 0), measured: { height: 100.4 } }]);
    expect(a).toBe(b);
  });
});
