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
