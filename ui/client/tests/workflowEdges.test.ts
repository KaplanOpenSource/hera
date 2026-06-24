import { describe, it, expect } from 'vitest';
import { buildWorkflowEdges, isValidConnection } from '../src/components/workflow/workflowEdges';

describe('buildWorkflowEdges', () => {
  it('returns no edges when nothing has requires', () => {
    expect(buildWorkflowEdges(['a', 'b'], { a: {}, b: {} })).toEqual([]);
  });

  it('builds an edge from a single requires', () => {
    expect(buildWorkflowEdges(['a', 'b'], { a: {}, b: { requires: 'a' } })).toEqual([
      { id: 'a->b', source: 'a', target: 'b' },
    ]);
  });

  it('builds one edge per entry in a requires array', () => {
    const edges = buildWorkflowEdges(['a', 'b', 'c'], { a: {}, b: {}, c: { requires: ['a', 'b'] } });
    expect(edges).toEqual([
      { id: 'a->c', source: 'a', target: 'c' },
      { id: 'b->c', source: 'b', target: 'c' },
    ]);
  });

  it('skips requires that point to unknown nodes', () => {
    expect(buildWorkflowEdges(['a'], { a: { requires: 'ghost' } })).toEqual([]);
  });
});

describe('isValidConnection', () => {
  const nodes = { a: {}, b: { requires: 'a' } };
  const names = ['a', 'b'];

  it('rejects a self connection', () => {
    expect(isValidConnection({ source: 'a', target: 'a' }, names, nodes)).toBe(false);
  });

  it('rejects a missing endpoint', () => {
    expect(isValidConnection({ source: '', target: 'b' }, names, nodes)).toBe(false);
  });

  it('rejects a duplicate edge', () => {
    // b already requires a.
    expect(isValidConnection({ source: 'a', target: 'b' }, names, nodes)).toBe(false);
  });

  it('rejects an edge that would create a cycle', () => {
    // a→b exists; making a require b closes a loop.
    expect(isValidConnection({ source: 'b', target: 'a' }, names, nodes)).toBe(false);
  });

  it('accepts a valid new edge', () => {
    const plain = { a: {}, b: {} };
    expect(isValidConnection({ source: 'a', target: 'b' }, ['a', 'b'], plain)).toBe(true);
  });
});
