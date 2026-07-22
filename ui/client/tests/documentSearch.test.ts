import { describe, it, expect } from 'vitest';
import {
  parseSearchQuery,
  unknownSearchFields,
  documentMatchesQuery,
  documentSearchText,
  SearchTermKind,
} from '../src/utils/documentSearch';

// Minimal document-data shape the search runs against.
const doc = {
  type: 'notebook',
  resource: '/projects/urban/flat.ipynb',
  desc: {
    groupName: 'flatUrban',
    datasourceName: 'analysis',
    effectParameters: { tenbergeCoefficient: 10 },
  },
};

const matches = (query: string): boolean => {
  const terms = parseSearchQuery(query);
  return documentMatchesQuery(doc, documentSearchText(doc), terms);
};

describe('parseSearchQuery', () => {
  it('classifies a bare word as a plain term', () => {
    expect(parseSearchQuery('flat')).toEqual([{ kind: SearchTermKind.Plain, value: 'flat' }]);
  });

  it('classifies a known root as a field term with a resolved path', () => {
    expect(parseSearchQuery('desc.groupName:flat')).toEqual([
      { kind: SearchTermKind.Field, value: 'flat', path: ['desc', 'groupName'] },
    ]);
  });

  it('lowercases and canonicalises the root but keeps sub-path casing', () => {
    expect(parseSearchQuery('Desc.groupName:Flat')).toEqual([
      { kind: SearchTermKind.Field, value: 'flat', path: ['desc', 'groupName'] },
    ]);
  });

  it('classifies an unknown root as an unknown term', () => {
    expect(parseSearchQuery('foo:bar')).toEqual([
      { kind: SearchTermKind.Unknown, value: 'bar', field: 'foo' },
    ]);
  });

  it('treats a value with a colon (non-identifier left side) as plain text', () => {
    expect(parseSearchQuery('12:30')).toEqual([{ kind: SearchTermKind.Plain, value: '12:30' }]);
  });

  it('splits on whitespace into multiple terms', () => {
    expect(parseSearchQuery('type:notebook flat')).toHaveLength(2);
  });
});

describe('unknownSearchFields', () => {
  it('collects and deduplicates unknown field roots', () => {
    expect(unknownSearchFields(parseSearchQuery('foo:1 foo:2 bar:3'))).toEqual(['foo', 'bar']);
  });

  it('is empty when all fields are known', () => {
    expect(unknownSearchFields(parseSearchQuery('type:notebook flat'))).toEqual([]);
  });
});

describe('documentMatchesQuery', () => {
  it('matches a plain term against any value (substring)', () => {
    expect(matches('urban')).toBe(true);
    expect(matches('missing')).toBe(false);
  });

  it('matches a scoped field term by substring', () => {
    expect(matches('type:note')).toBe(true);
    expect(matches('type:sim')).toBe(false);
  });

  it('drills into a dotted path', () => {
    expect(matches('desc.groupName:flat')).toBe(true);
    expect(matches('desc.effectParameters.tenbergeCoefficient:10')).toBe(true);
    expect(matches('desc.groupName:coast')).toBe(false);
  });

  it('ANDs multiple terms', () => {
    expect(matches('type:notebook desc.groupName:flat')).toBe(true);
    expect(matches('type:notebook desc.groupName:coast')).toBe(false);
  });

  it('never matches an unknown field', () => {
    expect(matches('foo:bar')).toBe(false);
  });

  it('matches everything for an empty query', () => {
    expect(matches('')).toBe(true);
  });
});
