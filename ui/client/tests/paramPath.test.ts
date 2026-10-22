import { describe, it, expect } from 'vitest';
import {
  leafStringsOf,
  paramPathOf,
  segmentsOf,
  valueAtPath,
  withoutPath,
  withValueAtPath,
} from '../src/components/workflow/paramPath';

describe('paramPathOf and segmentsOf', () => {
  it('joins and splits on dots', () => {
    expect(paramPathOf(['P', 'one'])).toBe('P.one');
    expect(segmentsOf('P.one')).toEqual(['P', 'one']);
  });

  it('round-trips a top-level name', () => {
    expect(segmentsOf(paramPathOf(['cmd']))).toEqual(['cmd']);
  });
});

describe('valueAtPath', () => {
  const params = { cmd: 'ls', P: { one: 'a', deep: { two: 'b' } }, list: ['x', 'y'] };

  it('reads a top-level value', () => {
    expect(valueAtPath(params, 'cmd')).toBe('ls');
  });

  it('reads a value inside a dict', () => {
    expect(valueAtPath(params, 'P.one')).toBe('a');
    expect(valueAtPath(params, 'P.deep.two')).toBe('b');
  });

  it('reads a list element by its position', () => {
    expect(valueAtPath(params, 'list.1')).toBe('y');
  });

  it('is undefined for a path that is not there', () => {
    expect(valueAtPath(params, 'nope')).toBeUndefined();
    expect(valueAtPath(params, 'P.nope.deeper')).toBeUndefined();
    expect(valueAtPath(params, 'cmd.one')).toBeUndefined();
  });
});

describe('withValueAtPath', () => {
  it('sets a top-level value', () => {
    expect(withValueAtPath({ a: '1' }, 'b', '2')).toEqual({ a: '1', b: '2' });
  });

  it('sets a value inside a dict and keeps its siblings', () => {
    expect(withValueAtPath({ P: { one: 'a', two: 'b' } }, 'P.one', 'z'))
      .toEqual({ P: { one: 'z', two: 'b' } });
  });

  it('sets a list element and keeps it a list', () => {
    const updated = withValueAtPath({ list: ['x', 'y'] }, 'list.1', 'z');
    expect(updated.list).toEqual(['x', 'z']);
    expect(Array.isArray(updated.list)).toBe(true);
  });

  it('creates the levels on the way when they are missing', () => {
    expect(withValueAtPath({}, 'P.one', 'z')).toEqual({ P: { one: 'z' } });
  });

  it('leaves the original untouched', () => {
    const params = { P: { one: 'a' } };
    withValueAtPath(params, 'P.one', 'z');
    expect(params.P.one).toBe('a');
  });
});

describe('withoutPath', () => {
  it('drops a top-level parameter', () => {
    expect(withoutPath({ a: '1', b: '2' }, 'a')).toEqual({ b: '2' });
  });

  it('drops a key inside a dict', () => {
    expect(withoutPath({ P: { one: 'a', two: 'b' } }, 'P.one')).toEqual({ P: { two: 'b' } });
  });

  it('splices a list element out', () => {
    expect(withoutPath({ list: ['x', 'y', 'z'] }, 'list.1')).toEqual({ list: ['x', 'z'] });
  });

  it('changes nothing for a path that is not there', () => {
    expect(withoutPath({ a: '1' }, 'P.one')).toEqual({ a: '1' });
  });

  it('leaves the original untouched', () => {
    const params = { P: { one: 'a' } };
    withoutPath(params, 'P.one');
    expect(params.P.one).toBe('a');
  });
});

describe('leafStringsOf', () => {
  it('returns a top-level string with its name as the path', () => {
    expect(leafStringsOf({ cmd: 'ls' })).toEqual([{ path: 'cmd', text: 'ls' }]);
  });

  it('walks into a dict', () => {
    expect(leafStringsOf({ P: { one: 'a', two: 'b' } })).toEqual([
      { path: 'P.one', text: 'a' },
      { path: 'P.two', text: 'b' },
    ]);
  });

  it('walks into a list, by position', () => {
    expect(leafStringsOf({ Command: ['a', 'b'] })).toEqual([
      { path: 'Command.0', text: 'a' },
      { path: 'Command.1', text: 'b' },
    ]);
  });

  it('walks a dict inside a list', () => {
    expect(leafStringsOf({ items: [{ k: 'v' }] })).toEqual([{ path: 'items.0.k', text: 'v' }]);
  });

  it('skips values that are not strings', () => {
    expect(leafStringsOf({ n: 5, b: true, z: null, s: 'yes' })).toEqual([{ path: 's', text: 'yes' }]);
  });

  it('is empty for no parameters', () => {
    expect(leafStringsOf({})).toEqual([]);
  });
});
