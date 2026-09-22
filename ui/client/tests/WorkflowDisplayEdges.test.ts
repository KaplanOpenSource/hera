import { describe, it, expect, vi } from 'vitest';
import { WorkflowDisplayEdges } from '../src/components/workflow/WorkflowDisplayEdges';

const requiresEdge = { id: 'A->B', source: 'A', target: 'B' };
const dataflowEdge = {
  id: 'df:A.out->B.param',
  source: 'A',
  sourceHandle: 'A:out:out',
  target: 'B',
  targetHandle: 'B:in:param',
};

describe('WorkflowDisplayEdges', () => {
  it('starts empty', () => {
    expect(WorkflowDisplayEdges.hovering(null).all()).toEqual([]);
  });

  it('attaches the requires handles', () => {
    const edges = WorkflowDisplayEdges.hovering(null).withRequires([requiresEdge], vi.fn()).all();
    expect(edges[0].type).toBe('requires');
    expect(edges[0].sourceHandle).toBe('A:req-out');
    expect(edges[0].targetHandle).toBe('B:req-in');
  });

  it('keeps both kinds, requires first', () => {
    const edges = WorkflowDisplayEdges.hovering(null)
      .withRequires([requiresEdge], vi.fn())
      .withDataflow([dataflowEdge], 'blue', vi.fn())
      .all();
    expect(edges.map((edge) => { return edge.type; })).toEqual(['requires', 'dataflow']);
  });

  it('marks only the hovered edge', () => {
    const edges = WorkflowDisplayEdges.hovering('A->B')
      .withRequires([requiresEdge], vi.fn())
      .withDataflow([dataflowEdge], 'blue', vi.fn())
      .all();
    expect(edges.map((edge) => { return edge.data!.hovered; })).toEqual([true, false]);
  });

  it('removes a requires edge by its ends', () => {
    const onRemove = vi.fn();
    const edges = WorkflowDisplayEdges.hovering(null).withRequires([requiresEdge], onRemove).all();
    (edges[0].data!.onRemove as () => void)();
    expect(onRemove).toHaveBeenCalledWith('A', 'B');
  });

  it('removes a dataflow edge by its id', () => {
    const onRemove = vi.fn();
    const edges = WorkflowDisplayEdges.hovering(null).withDataflow([dataflowEdge], 'blue', onRemove).all();
    (edges[0].data!.onRemove as () => void)();
    expect(onRemove).toHaveBeenCalledWith(dataflowEdge.id);
  });

  it('colors the dataflow line and its arrow', () => {
    const edges = WorkflowDisplayEdges.hovering(null).withDataflow([dataflowEdge], 'blue', vi.fn()).all();
    expect(edges[0].style!.stroke).toBe('blue');
    expect((edges[0].markerEnd as { color: string }).color).toBe('blue');
  });

  it('leaves the earlier instance untouched', () => {
    const base = WorkflowDisplayEdges.hovering(null).withRequires([requiresEdge], vi.fn());
    base.withDataflow([dataflowEdge], 'blue', vi.fn());
    expect(base.all()).toHaveLength(1);
  });
});
