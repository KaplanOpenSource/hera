import { describe, it, expect } from 'vitest';
import {
  nodeWithoutParam,
  nodeWithParamValue,
  nodeWithReferenceAt,
  nodeWithRequire,
  nodeWithRequiresFromText,
} from '../src/components/workflow/workflowNodeEdits';
import { Reference } from '../src/components/workflow/references/Reference';
import { OUTPUT } from '../src/components/workflow/references/knownKinds';
import { WorkflowNode } from '../src/shared/types';

const node = (params: { [param: string]: any }): WorkflowNode => {
  return { type: 'maker', Execution: { input_parameters: params } };
};

describe('nodeWithoutParam', () => {
  it('drops the named parameter and keeps the rest', () => {
    const updated = nodeWithoutParam(node({ a: '1', b: '2' }), 'a');
    expect(updated.Execution?.input_parameters).toEqual({ b: '2' });
  });

  it('is a no-op for a parameter that is not there', () => {
    const updated = nodeWithoutParam(node({ a: '1' }), 'zzz');
    expect(updated.Execution?.input_parameters).toEqual({ a: '1' });
  });

  it('leaves the original node untouched', () => {
    const original = node({ a: '1' });
    nodeWithoutParam(original, 'a');
    expect(original.Execution?.input_parameters).toEqual({ a: '1' });
  });

  it('copes with a node that has no Execution at all', () => {
    expect(nodeWithoutParam({ type: 'maker' }, 'a').Execution?.input_parameters).toEqual({});
  });

  it('drops a key inside a dict parameter', () => {
    const updated = nodeWithoutParam(node({ P: { one: '1', two: '2' } }), 'P.one');
    expect(updated.Execution?.input_parameters).toEqual({ P: { two: '2' } });
  });
});

describe('nodeWithParamValue', () => {
  it('sets a new parameter', () => {
    expect(nodeWithParamValue(node({}), 'a', '1').Execution?.input_parameters).toEqual({ a: '1' });
  });

  it('overwrites an existing one and keeps the others', () => {
    const updated = nodeWithParamValue(node({ a: '1', b: '2' }), 'a', '9');
    expect(updated.Execution?.input_parameters).toEqual({ a: '9', b: '2' });
  });

  it('keeps a non-string value as it is', () => {
    expect(nodeWithParamValue(node({}), 'a', 5).Execution?.input_parameters?.a).toBe(5);
  });

  it('keeps the rest of the node', () => {
    expect(nodeWithParamValue(node({}), 'a', '1').type).toBe('maker');
  });

  it('sets a value inside a dict parameter', () => {
    const updated = nodeWithParamValue(node({ P: { one: '1' } }), 'P.one', '9');
    expect(updated.Execution?.input_parameters).toEqual({ P: { one: '9' } });
  });
});

describe('nodeWithReferenceAt', () => {
  it('appends the reference when there is no caret', () => {
    const updated = nodeWithReferenceAt(node({ cmd: 'run' }), 'cmd', new Reference('A', OUTPUT, 'result'));
    expect(updated.Execution?.input_parameters?.cmd).toBe('run{A.output.result}');
  });

  it('inserts at the caret, keeping the text on both sides', () => {
    const updated = nodeWithReferenceAt(node({ cmd: 'ab' }), 'cmd', new Reference('A', OUTPUT, 'result'), 1);
    expect(updated.Execution?.input_parameters?.cmd).toBe('a{A.output.result}b');
  });

  it('starts from empty text when the parameter is missing', () => {
    const updated = nodeWithReferenceAt(node({}), 'cmd', new Reference('A', OUTPUT, 'result'));
    expect(updated.Execution?.input_parameters?.cmd).toBe('{A.output.result}');
  });

  it('replaces a non-string value rather than appending to it', () => {
    const updated = nodeWithReferenceAt(node({ cmd: 5 }), 'cmd', new Reference('A', OUTPUT, 'result'));
    expect(updated.Execution?.input_parameters?.cmd).toBe('{A.output.result}');
  });

  it('inserts into a value inside a dict parameter', () => {
    const updated = nodeWithReferenceAt(node({ P: { one: 'ab' } }), 'P.one', new Reference('A', OUTPUT, 'result'), 1);
    expect(updated.Execution?.input_parameters?.P).toEqual({ one: 'a{A.output.result}b' });
  });

  it('leaves the other parameters alone', () => {
    const updated = nodeWithReferenceAt(node({ cmd: '', keep: 'as is' }), 'cmd', new Reference('A', OUTPUT, 'result'));
    expect(updated.Execution?.input_parameters?.keep).toBe('as is');
  });
});

describe('nodeWithRequire', () => {
  it('adds the name to requires', () => {
    expect(nodeWithRequire(node({}), 'A', 'B').requires).toEqual(['A']);
  });

  it('keeps a name that is already required', () => {
    const before = { ...node({}), requires: ['A'] };
    expect(nodeWithRequire(before, 'A', 'B').requires).toEqual(['A']);
  });

  it('turns a single name into a list', () => {
    const before = { ...node({}), requires: 'A' };
    expect(nodeWithRequire(before, 'C', 'B').requires).toEqual(['A', 'C']);
  });

  it('never makes a node require itself', () => {
    expect(nodeWithRequire(node({}), 'B', 'B').requires).toBeUndefined();
  });
});

describe('nodeWithRequiresFromText', () => {
  it('requires every node the text references', () => {
    const updated = nodeWithRequiresFromText(node({}), 'C', '{A.output.x} {B.output.y}', ['A', 'B', 'C']);
    expect(updated.requires).toEqual(['A', 'B']);
  });

  it('ignores a name that is not a node of this workflow', () => {
    const updated = nodeWithRequiresFromText(node({}), 'C', '{Zz.output.x}', ['A', 'C']);
    expect(updated.requires).toBeUndefined();
  });

  it('ignores text with no reference in it', () => {
    expect(nodeWithRequiresFromText(node({}), 'C', 'plain text', ['A', 'C']).requires).toBeUndefined();
  });
});
