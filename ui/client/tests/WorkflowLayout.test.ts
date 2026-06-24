import { describe, it, expect } from 'vitest';
import { WorkflowLayout } from '../src/components/workflow/WorkflowLayout';
import { V_GAP, X_GAP, estimateHeight } from '../src/components/workflow/workflowGeometry';

describe('WorkflowLayout.stacked', () => {
  it('stacks independent nodes in one column', () => {
    const layout = WorkflowLayout.stacked(['a', 'b'], { a: {}, b: {} }).positions();
    expect(layout.a).toEqual({ x: 0, y: 0 });
    expect(layout.b).toEqual({ x: 0, y: estimateHeight({}) + V_GAP });
  });

  it('places dependent nodes in successive columns', () => {
    const layout = WorkflowLayout.stacked(['a', 'b'], { a: {}, b: { requires: 'a' } }).positions();
    expect(layout.a).toEqual({ x: 0, y: 0 });
    expect(layout.b).toEqual({ x: X_GAP, y: 0 });
  });
});

describe('WorkflowLayout.fixOverlaps', () => {
  const fix = (placed: { id: string, layer: number, x: number, y: number, height: number }[], vGap: number) => {
    return WorkflowLayout.fromPlaced(placed).fixOverlaps(vGap).positions();
  };

  it('leaves non-overlapping nodes untouched', () => {
    const layout = fix([
      { id: 'a', layer: 0, x: 0, y: 0, height: 100 },
      { id: 'b', layer: 0, x: 0, y: 500, height: 100 },
    ], 20);
    expect(layout.a.y).toBe(0);
    expect(layout.b.y).toBe(500);
  });

  it('pushes a colliding node down to clear the one above plus the gap', () => {
    const layout = fix([
      { id: 'a', layer: 0, x: 0, y: 0, height: 200 },
      { id: 'b', layer: 0, x: 0, y: 50, height: 100 },
    ], 20);
    // a occupies 0..200; b must start at 200 + gap(20) = 220.
    expect(layout.a.y).toBe(0);
    expect(layout.b.y).toBe(220);
  });

  it('never moves a node up and never moves the topmost node', () => {
    const layout = fix([
      { id: 'a', layer: 0, x: 0, y: 100, height: 50 },
      { id: 'b', layer: 0, x: 0, y: 400, height: 50 },
    ], 20);
    expect(layout.a.y).toBe(100);
    expect(layout.b.y).toBe(400);
  });

  it('cascades a push through several stacked nodes', () => {
    const layout = fix([
      { id: 'a', layer: 0, x: 0, y: 0, height: 300 },
      { id: 'b', layer: 0, x: 0, y: 100, height: 100 },
      { id: 'c', layer: 0, x: 0, y: 250, height: 100 },
    ], 10);
    // a: 0..300 → b to 310 (310..410) → c to b.bottom 420.
    expect(layout.b.y).toBe(310);
    expect(layout.c.y).toBe(420);
  });

  it('resolves each column independently and keeps x', () => {
    const layout = fix([
      { id: 'a', layer: 0, x: 0, y: 0, height: 200 },
      { id: 'b', layer: 0, x: 0, y: 50, height: 100 },
      { id: 'c', layer: 1, x: X_GAP, y: 0, height: 100 },
    ], 20);
    // Column 0 collides (b pushed); column 1 has a single node, untouched.
    expect(layout.b).toEqual({ x: 0, y: 220 });
    expect(layout.c).toEqual({ x: X_GAP, y: 0 });
  });
});
