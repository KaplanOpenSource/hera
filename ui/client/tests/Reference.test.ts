import { describe, expect, it } from 'vitest';
import { Reference } from '../src/components/workflow/references/Reference';
import { OUTPUT } from '../src/components/workflow/references/knownKinds';

const reference = new Reference('A', OUTPUT, 'result');

describe('Reference', () => {
  it('writes the token a parameter value holds', () => {
    expect(reference.toString()).toBe('{A.output.result}');
  });

  it('names the source dot a line leaves from', () => {
    expect(reference.handleId()).toBe('A:out:result');
  });

  it('names the canvas line to the parameter reading it', () => {
    expect(reference.edgeIdTo('B', 'cmd')).toBe('df:A:out:result->B.cmd');
  });

  it('strips its own token and leaves the others', () => {
    const value = '{A.output.result} {A.output.log}';
    expect(value.replace(reference.clearToken(), '').trim()).toBe('{A.output.log}');
  });
});

describe('a key that points inside the output', () => {
  const deep = new Reference('A', OUTPUT, 'items[0].name');

  it('writes the whole path into the token', () => {
    expect(deep.toString()).toBe('{A.output.items[0].name}');
  });

  it('names the output itself, without the path', () => {
    expect(deep.rootKey()).toBe('items');
    expect(new Reference('A', OUTPUT, 'result.station').rootKey()).toBe('result');
    expect(reference.rootKey()).toBe('result');
  });

  it('leaves from the dot of the output it points into', () => {
    expect(deep.rootHandleId()).toBe('A:out:items');
  });

  it('strips its own token and leaves a sibling path alone', () => {
    const value = '{A.output.items[0].name} {A.output.items[1].name}';
    expect(value.replace(deep.clearToken(), '').trim()).toBe('{A.output.items[1].name}');
  });
});
