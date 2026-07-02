import { describe, it, expect } from 'vitest';
import { NodeCatalogEntry, nodeTypeGroup, prefilledParameters, nodeTypeIssue, paramsFieldDef, nodeOutputNames } from '../src/components/workflow/nodeCatalog';

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

describe('nodeTypeIssue', () => {
  const catalog = [entry('openFOAM.constant.g', ['x'])];

  it('reports nothing while the catalog is empty (still loading)', () => {
    expect(nodeTypeIssue({ type: 'anything' }, [])).toBeUndefined();
  });

  it('flags a node with no type', () => {
    expect(nodeTypeIssue({ type: '' }, catalog)).toBe('no type selected');
  });

  it('flags an unknown type', () => {
    expect(nodeTypeIssue({ type: 'made.up' }, catalog)).toBe('unknown type: made.up');
  });

  it('is fine for a known type', () => {
    expect(nodeTypeIssue({ type: 'openFOAM.constant.g' }, catalog)).toBeUndefined();
  });
});

describe('paramsFieldDef', () => {
  const catalog = [entry('openFOAM.constant.g', ['x', 'y'])];

  it('puts each parameter under children, marked required or not', () => {
    expect(paramsFieldDef({ type: 'openFOAM.constant.g' }, catalog)).toEqual({
      children: { x: { required: true }, y: { required: true } },
    });
  });

  it('reflects optional parameters as not required', () => {
    const cat: NodeCatalogEntry[] = [{
      type: 't',
      parameters: [
        { name: 'req', is_required: true, source: 'jsonForm' },
        { name: 'opt', is_required: false, source: 'template' },
      ],
    }];
    expect(paramsFieldDef({ type: 't' }, cat)).toEqual({
      children: { req: { required: true }, opt: { required: false } },
    });
  });

  it('has empty children for an unknown type', () => {
    expect(paramsFieldDef({ type: 'made.up' }, catalog)).toEqual({ children: {} });
  });
});

describe('nodeOutputNames', () => {
  const catalog: NodeCatalogEntry[] = [{
    type: 'general.CopyFile',
    parameters: [],
    outputs: [
      { name: 'case', source: 'python' },
      { name: 'fields', source: 'python' },
    ],
  }];

  it('returns the output names for a known type', () => {
    expect(nodeOutputNames({ type: 'general.CopyFile' }, catalog)).toEqual(['case', 'fields']);
  });

  it('returns empty for an unknown type', () => {
    expect(nodeOutputNames({ type: 'made.up' }, catalog)).toEqual([]);
  });

  it('returns empty when the type has no outputs field', () => {
    expect(nodeOutputNames({ type: 't' }, [{ type: 't', parameters: [] }])).toEqual([]);
  });
});
