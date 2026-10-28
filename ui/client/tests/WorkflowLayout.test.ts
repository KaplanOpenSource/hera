import { describe, it, expect } from 'vitest';
import { WorkflowLayout } from '../src/components/workflow/WorkflowLayout';
import { H_GAP, V_GAP, X_GAP, estimateHeight, estimatedWidth } from '../src/components/workflow/workflowGeometry';

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
    return WorkflowLayout.fromPlaced(placed.map(node => ({ ...node, width: 300 }))).fixOverlaps(vGap).positions();
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

describe('WorkflowLayout.compact', () => {
  // Every node 300 wide unless the test says otherwise, so the column x values
  // are easy to read: 0, then 300 + hGap, ...
  const compact = (
    placed: { id: string, layer: number, x: number, y: number, height: number, width?: number }[],
    vGap: number,
    hGap: number = H_GAP,
  ) => {
    return WorkflowLayout.fromPlaced(placed.map(node => ({ width: 300, ...node }))).compact(vGap, hGap).positions();
  };

  it('pulls a node up to one gap under the one above', () => {
    const layout = compact([
      { id: 'a', layer: 0, x: 0, y: 0, height: 100 },
      { id: 'b', layer: 0, x: 0, y: 500, height: 100 },
    ], 20);
    expect(layout.a.y).toBe(0);
    expect(layout.b.y).toBe(120);
  });

  it('pushes a node down when the one above grew', () => {
    const layout = compact([
      { id: 'a', layer: 0, x: 0, y: 0, height: 200 },
      { id: 'b', layer: 0, x: 0, y: 50, height: 100 },
    ], 20);
    expect(layout.b.y).toBe(220);
  });

  it('leaves the topmost node of a column where it is', () => {
    const layout = compact([
      { id: 'a', layer: 0, x: 0, y: 300, height: 100 },
      { id: 'b', layer: 0, x: 0, y: 900, height: 100 },
    ], 20);
    expect(layout.a.y).toBe(300);
    expect(layout.b.y).toBe(420);
  });

  it('keeps the order the nodes are already in', () => {
    const layout = compact([
      { id: 'a', layer: 0, x: 0, y: 900, height: 100 },
      { id: 'b', layer: 0, x: 0, y: 0, height: 100 },
    ], 20);
    expect(layout.b.y).toBe(0);
    expect(layout.a.y).toBe(120);
  });

  it('compacts each column on its own', () => {
    const layout = compact([
      { id: 'a', layer: 0, x: 0, y: 0, height: 100 },
      { id: 'b', layer: 0, x: 0, y: 800, height: 100 },
      { id: 'c', layer: 1, x: X_GAP, y: 0, height: 100 },
    ], 20, 40);
    expect(layout.b).toEqual({ x: 0, y: 120 });
    expect(layout.c).toEqual({ x: 340, y: 0 });
  });

  it('puts a column one gap right of the widest node before it', () => {
    const layout = compact([
      { id: 'a', layer: 0, x: 0, y: 0, height: 100, width: 300 },
      { id: 'b', layer: 0, x: 0, y: 200, height: 100, width: 560 },
      { id: 'c', layer: 1, x: 9999, y: 0, height: 100, width: 300 },
    ], 20, 40);
    // The widest node of column 0 is 560, so column 1 starts at 600.
    expect(layout.c.x).toBe(600);
  });

  it('pulls a column back left when the node before it shrank', () => {
    const layout = compact([
      { id: 'a', layer: 0, x: 0, y: 0, height: 100, width: 300 },
      { id: 'b', layer: 1, x: 5000, y: 0, height: 100, width: 300 },
    ], 20, 40);
    expect(layout.b.x).toBe(340);
  });

  it('keeps the leftmost column where it is', () => {
    const layout = compact([
      { id: 'a', layer: 0, x: 700, y: 0, height: 100, width: 300 },
      { id: 'b', layer: 1, x: 0, y: 0, height: 100, width: 300 },
    ], 20, 40);
    expect(layout.a.x).toBe(700);
    expect(layout.b.x).toBe(1040);
  });
});

describe('estimatedWidth', () => {
  it('leaves the same room as the estimated column pitch', () => {
    expect(estimatedWidth() + H_GAP).toBe(X_GAP);
  });
});
