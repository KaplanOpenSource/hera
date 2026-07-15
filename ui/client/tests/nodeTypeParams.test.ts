import { describe, it, expect } from 'vitest';
import { NodeCatalogEntry } from '../src/components/workflow/nodeCatalog';
import { paramsOnTypeChange } from '../src/components/workflow/nodeTypeParams';
import { NodeParameterSource } from '../src/shared/types';

const entry = (type: string, names: string[]): NodeCatalogEntry => ({
  type,
  parameters: names.map(name => ({ name, is_required: true, source: NodeParameterSource.JsonForm })),
});

describe('paramsOnTypeChange', () => {
  it('seeds the new type params as empty when nothing exists', () => {
    expect(paramsOnTypeChange({}, entry('t', ['transportModel', 'nu']))).toEqual({ transportModel: '', nu: '' });
  });

  it('keeps a param that has a value and adds the missing ones', () => {
    expect(paramsOnTypeChange({ a: 5 }, entry('t', ['a', 'b']))).toEqual({ a: 5, b: '' });
  });

  it('drops empty params of the old type and seeds the new type', () => {
    expect(paramsOnTypeChange({ old: '', keep: 'hi' }, entry('t', ['x', 'y'])))
      .toEqual({ keep: 'hi', x: '', y: '' });
  });

  it('keeps a shared param value (e.g. ProjectName) across the change', () => {
    expect(paramsOnTypeChange({ ProjectName: 'demo' }, entry('t', ['ProjectName', 'z'])))
      .toEqual({ ProjectName: 'demo', z: '' });
  });

  it('keeps 0 and false — they are real values, not empty', () => {
    expect(paramsOnTypeChange({ n: 0, b: false }, entry('t', ['x']))).toEqual({ n: 0, b: false, x: '' });
  });

  it('skips nested/template paths and loop variables', () => {
    expect(paramsOnTypeChange({}, entry('t', ['fields', 'fields.items', 'a b']))).toEqual({ fields: '' });
  });

  it('returns the existing params unchanged for an unknown type', () => {
    const existing = { a: 1 };
    expect(paramsOnTypeChange(existing, undefined)).toBe(existing);
  });
});
