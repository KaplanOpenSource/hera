import { describe, expect, it } from 'vitest';
import { InputReferenceKind } from '../src/components/workflow/references/InputReferenceKind';
import { Reference } from '../src/components/workflow/references/Reference';

const INPUT = new InputReferenceKind();

describe('InputReferenceKind', () => {
  it('writes the token a parameter value holds', () => {
    expect(new Reference('A', INPUT, 'cmd').toString()).toBe('{A.Execution.input_parameters.cmd}');
  });

  it('names the source dot a line leaves from', () => {
    expect(new Reference('A', INPUT, 'cmd').handleId()).toBe('A:param:cmd');
  });

  it('lists the node own parameter keys', () => {
    const node = { Execution: { input_parameters: { cmd: 'ls', other: 1 } } };
    expect(INPUT.namesOf(node, [])).toEqual(['cmd', 'other']);
  });

  it('lists nothing for a node with no parameters', () => {
    expect(INPUT.namesOf({}, [])).toEqual([]);
  });
});
