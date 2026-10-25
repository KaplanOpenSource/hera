import { describe, it, expect } from 'vitest';
import {
  BASE_HEIGHT,
  ROW_HEIGHT,
  computeLayers,
  estimateHeight,
} from '../src/components/workflow/workflowGeometry';

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

describe('estimateHeight', () => {
  it('is just the base for a node with no params', () => {
    expect(estimateHeight({})).toBe(BASE_HEIGHT);
  });

  it('grows by a row per parameter', () => {
    const node = { Execution: { input_parameters: { a: 1, b: 2 } } };
    expect(estimateHeight(node)).toBe(BASE_HEIGHT + 2 * ROW_HEIGHT);
  });

  // Nested values are rows of the editor, not of the summary card.
  it('counts a nested dict as one row', () => {
    const node = { Execution: { input_parameters: { a: { b: 1, c: 2 } } } };
    expect(estimateHeight(node)).toBe(BASE_HEIGHT + ROW_HEIGHT);
  });
});

