import { describe, expect, it } from 'vitest';
import { NodeParameterSource } from '../src/shared/types';
import { NodeCatalogEntry } from '../src/components/workflow/nodeCatalog';
import { Reference } from '../src/components/workflow/references/Reference';
import { INPUT, OUTPUT, knownKinds } from '../src/components/workflow/references/knownKinds';

const catalog: NodeCatalogEntry[] = [{
  type: 'general.Run',
  parameters: [],
  outputs: [
    { name: 'result', source: NodeParameterSource.Python },
    { name: 'log', source: NodeParameterSource.Python },
  ],
}];

describe('parseAll', () => {
  it('parses back what a reference wrote', () => {
    const reference = new Reference('A', OUTPUT, 'result');
    expect(knownKinds.parseAll(reference.toString())).toEqual([reference]);
  });

  it('finds two references in one value', () => {
    expect(knownKinds.parseAll('run {A.output.result} then {C.output.log}')).toEqual([
      new Reference('A', OUTPUT, 'result'),
      new Reference('C', OUTPUT, 'log'),
    ]);
  });

  it('still reads the older parameters spelling', () => {
    expect(knownKinds.parseAll('{A.parameters.result}')).toEqual([new Reference('A', OUTPUT, 'result')]);
  });

  it('skips a section no kind knows', () => {
    expect(knownKinds.parseAll('{A.mystery.result}')).toEqual([]);
  });

  it('finds nothing in plain text', () => {
    expect(knownKinds.parseAll('just a command')).toEqual([]);
  });
});

describe('ofHandle', () => {
  it('round trips a source dot', () => {
    const reference = new Reference('A', OUTPUT, 'result');
    expect(knownKinds.ofHandle(reference.handleId())).toEqual(reference);
  });

  it('ignores a requires handle', () => {
    expect(knownKinds.ofHandle('A:req-out')).toBeNull();
  });

  it('ignores an input handle', () => {
    expect(knownKinds.ofHandle('A:in:cmd')).toBeNull();
  });
});

describe('ofEdgeId', () => {
  it('round trips a dataflow edge id', () => {
    const reference = new Reference('A', OUTPUT, 'result');
    const id = reference.edgeIdTo('B', 'cmd');
    expect(knownKinds.ofEdgeId(id)).toEqual({ reference, target: 'B', paramPath: 'cmd' });
  });

  it('ignores an id that is not a dataflow edge', () => {
    expect(knownKinds.ofEdgeId('A->B')).toBeNull();
  });
});

describe('the kinds themselves', () => {
  it('names the kind a written section belongs to', () => {
    expect(knownKinds.bySection('output')).toBe(OUTPUT);
    expect(knownKinds.bySection('outputs')).toBe(OUTPUT);
    expect(knownKinds.bySection('Execution.input_parameters')).toBe(INPUT);
    expect(knownKinds.bySection('mystery')).toBeNull();
  });

  it('keeps a half-typed section open', () => {
    expect(knownKinds.couldBe('out')).toEqual([OUTPUT]);
    expect(knownKinds.couldBe('Exec')).toEqual([INPUT]);
    expect(knownKinds.couldBe('zz')).toEqual([]);
  });

  it('lists what a node offers, of both kinds', () => {
    const node = { type: 'general.Run', Execution: { input_parameters: { cmd: 'ls' } } };
    expect(knownKinds.keysOf('A', node, catalog)).toEqual([
      new Reference('A', OUTPUT, 'result'),
      new Reference('A', OUTPUT, 'log'),
      new Reference('A', INPUT, 'cmd'),
    ]);
  });

  it('offers nothing for an unknown type', () => {
    expect(knownKinds.keysOf('A', { type: 'nope' }, catalog)).toEqual([]);
  });
});

describe('renamedNode', () => {
  it('renames an output reference', () => {
    expect(knownKinds.renamedNode('{a.output.x}', 'a', 'b')).toEqual('{b.output.x}');
  });

  it('renames an input reference', () => {
    expect(knownKinds.renamedNode('{a.Execution.input_parameters.k}', 'a', 'b'))
      .toEqual('{b.Execution.input_parameters.k}');
  });

  it('keeps the older parameters spelling', () => {
    expect(knownKinds.renamedNode('{a.parameters.x}', 'a', 'b')).toEqual('{b.parameters.x}');
  });

  it('leaves a different node alone', () => {
    expect(knownKinds.renamedNode('{c.output.x}', 'a', 'b')).toEqual('{c.output.x}');
  });

  it('leaves a name that only shares a prefix alone', () => {
    expect(knownKinds.renamedNode('{ab.output.x}', 'a', 'b')).toEqual('{ab.output.x}');
  });

  it('renames two tokens in one string', () => {
    expect(knownKinds.renamedNode('{a.output.x} {a.output.y}', 'a', 'b')).toEqual('{b.output.x} {b.output.y}');
  });

  it('leaves an unknown section alone', () => {
    expect(knownKinds.renamedNode('{a.nope.x}', 'a', 'b')).toEqual('{a.nope.x}');
  });

  it('drops whitespace inside the braces', () => {
    expect(knownKinds.renamedNode('{ a.output.x }', 'a', 'b')).toEqual('{b.output.x}');
  });

  it('keeps the text around the token', () => {
    expect(knownKinds.renamedNode('-{a.output.x}+1e-6', 'a', 'b')).toEqual('-{b.output.x}+1e-6');
  });
});

// Issue #1064: the key may be a JSONPath into the output, which Hermes resolves.
describe('a key that points inside the output', () => {
  it('parses a dict field', () => {
    expect(knownKinds.parseAll('{A.output.result.station}'))
      .toEqual([new Reference('A', OUTPUT, 'result.station')]);
  });

  it('parses a list element', () => {
    expect(knownKinds.parseAll('{A.output.items[0]}'))
      .toEqual([new Reference('A', OUTPUT, 'items[0]')]);
  });

  it('parses a field inside a list element', () => {
    expect(knownKinds.parseAll('{A.output.items[0].name}'))
      .toEqual([new Reference('A', OUTPUT, 'items[0].name')]);
  });

  it('parses every element of a list', () => {
    expect(knownKinds.parseAll('{A.output.items[*].name}'))
      .toEqual([new Reference('A', OUTPUT, 'items[*].name')]);
  });

  it('still reads an input reference, whose section holds a dot', () => {
    expect(knownKinds.parseAll('{A.Execution.input_parameters.cmd}'))
      .toEqual([new Reference('A', INPUT, 'cmd')]);
  });

  it('reads a sub-path on an input reference too', () => {
    expect(knownKinds.parseAll('{A.Execution.input_parameters.cmd.deep}'))
      .toEqual([new Reference('A', INPUT, 'cmd.deep')]);
  });

  it('skips a section with no key after it', () => {
    expect(knownKinds.parseAll('{A.output.}')).toEqual([]);
  });

  it('leaves braces that are not a reference alone', () => {
    expect(knownKinds.parseAll('echo {} and {a b}')).toEqual([]);
    expect(knownKinds.renamedNode('echo {} and {a b}', 'a', 'b')).toEqual('echo {} and {a b}');
  });

  it('renames a node on a reference with a sub-path', () => {
    expect(knownKinds.renamedNode('{a.output.items[0].name}', 'a', 'b'))
      .toEqual('{b.output.items[0].name}');
  });

  it('round trips an edge id whose key is a path', () => {
    const reference = new Reference('A', OUTPUT, 'items[0].name');
    expect(knownKinds.ofEdgeId(reference.edgeIdTo('B', 'Parameters.cmd')))
      .toEqual({ reference, target: 'B', paramPath: 'Parameters.cmd' });
  });

  it('round trips a source dot whose key is a path', () => {
    const reference = new Reference('A', OUTPUT, 'items[0].name');
    expect(knownKinds.ofHandle(reference.handleId())).toEqual(reference);
  });
});

describe('splitSection', () => {
  it('tells the section from the key', () => {
    expect(knownKinds.splitSection('output.result.station')).toEqual({ kind: OUTPUT, key: 'result.station' });
  });

  it('takes the longest section that matches', () => {
    expect(knownKinds.splitSection('Execution.input_parameters.cmd')).toEqual({ kind: INPUT, key: 'cmd' });
  });

  it('allows an empty key, for a section still being typed', () => {
    expect(knownKinds.splitSection('output.')).toEqual({ kind: OUTPUT, key: '' });
  });

  it('knows nothing of an unknown section', () => {
    expect(knownKinds.splitSection('mystery.x')).toBeNull();
  });
});
