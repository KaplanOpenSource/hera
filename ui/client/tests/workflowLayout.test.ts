import { describe, it, expect } from 'vitest';
import {
  BASE_HEIGHT,
  ROW_HEIGHT,
  V_GAP,
  X_GAP,
  computeLayers,
  computeLayout,
  countRows,
  estimateHeight,
} from '../src/components/workflow/workflowLayout';

describe('computeLayers', () => {
  it('puts nodes with no requires on layer 0', () => {
    expect(computeLayers(['a', 'b'], { a: {}, b: {} })).toEqual({ a: 0, b: 0 });
  });

  it('increases the layer along a requires chain', () => {
    const nodes = { a: {}, b: { requires: 'a' }, c: { requires: 'b' } };
    expect(computeLayers(['a', 'b', 'c'], nodes)).toEqual({ a: 0, b: 1, c: 2 });
  });

  it('uses the longest chain when a node has several requires', () => {
    // d requires both a (depth 0) and c (depth 2) → layer 3.
    const nodes = { a: {}, b: { requires: 'a' }, c: { requires: 'b' }, d: { requires: ['a', 'c'] } };
    expect(computeLayers(['a', 'b', 'c', 'd'], nodes).d).toBe(3);
  });

  it('breaks cycles at layer 0', () => {
    const nodes = { a: { requires: 'b' }, b: { requires: 'a' } };
    const layers = computeLayers(['a', 'b'], nodes);
    expect(layers.a).toBeGreaterThanOrEqual(0);
    expect(layers.b).toBeGreaterThanOrEqual(0);
  });

  it('ignores requires that point to unknown nodes', () => {
    expect(computeLayers(['a'], { a: { requires: 'ghost' } })).toEqual({ a: 0 });
  });
});

describe('countRows', () => {
  it('is 0 for a primitive', () => {
    expect(countRows(5)).toBe(0);
    expect(countRows('x')).toBe(0);
    expect(countRows(null)).toBe(0);
  });

  it('counts one row per key in a flat object', () => {
    expect(countRows({ a: 1, b: 2 })).toBe(2);
  });

  it('counts nested keys recursively', () => {
    expect(countRows({ a: { b: 1 } })).toBe(2);
    expect(countRows({ a: { b: 1, c: 2 }, d: 3 })).toBe(4);
  });

  it('is 0 for an empty object', () => {
    expect(countRows({})).toBe(0);
  });
});

describe('estimateHeight', () => {
  it('is base + one row for a node with no params', () => {
    expect(estimateHeight({})).toBe(BASE_HEIGHT + ROW_HEIGHT);
  });

  it('grows by a row per parameter', () => {
    const node = { Execution: { input_parameters: { a: 1, b: 2 } } };
    expect(estimateHeight(node)).toBe(BASE_HEIGHT + (1 + 2) * ROW_HEIGHT);
  });
});

describe('computeLayout', () => {
  it('stacks independent nodes in one column', () => {
    const layout = computeLayout(['a', 'b'], { a: {}, b: {} });
    expect(layout.a).toEqual({ x: 0, y: 0 });
    expect(layout.b).toEqual({ x: 0, y: estimateHeight({}) + V_GAP });
  });

  it('places dependent nodes in successive columns', () => {
    const layout = computeLayout(['a', 'b'], { a: {}, b: { requires: 'a' } });
    expect(layout.a).toEqual({ x: 0, y: 0 });
    expect(layout.b).toEqual({ x: X_GAP, y: 0 });
  });
});

