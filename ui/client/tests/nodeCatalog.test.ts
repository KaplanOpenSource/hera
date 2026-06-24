import { describe, it, expect } from 'vitest';
import { NodeCatalogEntry, nodeTypeGroup, prefilledParameters, validateNode } from '../src/components/workflow/nodeCatalog';

describe('nodeTypeGroup', () => {
  it('drops the last dotted segment', () => {
    expect(nodeTypeGroup('openFOAM.constant.physicalProperties')).toBe('openFOAM.constant');
    expect(nodeTypeGroup('openFOAM.mesh.BlockMesh')).toBe('openFOAM.mesh');
  });

  it('returns the whole type when it has no prefix', () => {
    expect(nodeTypeGroup('BC')).toBe('BC');
  });
});

const entry = (type: string, names: string[]): NodeCatalogEntry => ({
  type,
  parameters: names.map(name => ({ name, is_required: true, source: 'jsonForm' })),
});

describe('prefilledParameters', () => {
  it('seeds each parameter name with an empty value', () => {
    const result = prefilledParameters(entry('t', ['transportModel', 'nu']), {});
    expect(result).toEqual({ transportModel: '', nu: '' });
  });

  it('keeps existing values and only adds missing names', () => {
    const result = prefilledParameters(entry('t', ['a', 'b']), { a: 5 });
    expect(result).toEqual({ a: 5, b: '' });
  });

  it('skips nested/template paths and loop variables', () => {
    const result = prefilledParameters(entry('t', ['fields', 'fields.items', 'a b']), {});
    expect(result).toEqual({ fields: '' });
  });

  it('returns the existing params unchanged for an unknown type', () => {
    const existing = { a: 1 };
    expect(prefilledParameters(undefined, existing)).toBe(existing);
  });
});

describe('validateNode', () => {
  const catalog = [entry('openFOAM.constant.g', ['x', 'y', 'z'])];

  it('reports nothing while the catalog is empty (still loading)', () => {
    expect(validateNode({ type: 'anything' }, [])).toEqual([]);
  });

  it('flags a node with no type', () => {
    expect(validateNode({ type: '' }, catalog)).toEqual(['no type selected']);
  });

  it('flags an unknown type', () => {
    expect(validateNode({ type: 'made.up' }, catalog)).toEqual(['unknown type: made.up']);
  });

  it('flags each empty required parameter', () => {
    const node = { type: 'openFOAM.constant.g', Execution: { input_parameters: { x: 1, y: '' } } };
    expect(validateNode(node, catalog)).toEqual(['missing required: y', 'missing required: z']);
  });

  it('passes when all required parameters have values', () => {
    const node = { type: 'openFOAM.constant.g', Execution: { input_parameters: { x: 1, y: 2, z: 0 } } };
    expect(validateNode(node, catalog)).toEqual([]);
  });

  it('ignores optional parameters', () => {
    const cat: NodeCatalogEntry[] = [{
      type: 't',
      parameters: [{ name: 'opt', is_required: false, source: 'template' }],
    }];
    expect(validateNode({ type: 't', Execution: { input_parameters: {} } }, cat)).toEqual([]);
  });
});
