import { describe, expect, it } from 'vitest';
import { mapStrings } from '../src/utils/mapStrings';

const upper = (text: string) => {
  return text.toUpperCase();
};

describe('mapStrings', () => {
  it('rewrites a bare string', () => {
    expect(mapStrings('a', upper)).toEqual('A');
  });

  it('rewrites strings deep in dicts and lists', () => {
    expect(mapStrings({ a: { b: ['c', { d: 'e' }] } }, upper)).toEqual({ a: { b: ['C', { d: 'E' }] } });
  });

  it('leaves keys untouched', () => {
    expect(mapStrings({ a: 'a' }, upper)).toEqual({ a: 'A' });
  });

  it('keeps numbers, booleans and null as they were', () => {
    expect(mapStrings({ n: 1, b: true, z: null }, upper)).toEqual({ n: 1, b: true, z: null });
  });

  it('does not mutate the input', () => {
    const input = { a: { b: 'c' } };
    const copy = JSON.parse(JSON.stringify(input));
    mapStrings(input, upper);
    expect(input).toEqual(copy);
  });
});
