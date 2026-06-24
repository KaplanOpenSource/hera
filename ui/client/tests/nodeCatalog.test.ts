import { describe, it, expect } from 'vitest';
import { NodeCatalogEntry, nodeTypeGroup, prefilledParameters } from '../src/components/workflow/nodeCatalog';

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
