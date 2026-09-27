import { describe, it, expect } from 'vitest';
import { WorkflowReferences } from '../src/components/workflow/WorkflowReferences';
import { Reference } from '../src/components/workflow/references/Reference';
import { OUTPUT } from '../src/components/workflow/references/knownKinds';
import { NodeCatalogEntry } from '../src/components/workflow/nodeCatalog';
import { NodeParameterSource } from '../src/shared/types';

const output = (name: string) => {
  return { name, source: NodeParameterSource.JsonForm };
};

const catalog: NodeCatalogEntry[] = [
  { type: 'maker', parameters: [], outputs: [output('first'), output('second')] },
  { type: 'plain', parameters: [] },
];

const nodes = {
  A: { type: 'maker' },
  B: { type: 'plain' },
  C: { type: 'maker' },
};

const references = new WorkflowReferences(['A', 'B', 'C'], nodes, catalog);

// An output reference, the only kind there is today.
const outputRef = (node: string, key: string): Reference => new Reference(node, OUTPUT, key);

describe('WorkflowReferences', () => {
  it('lists the other nodes that offer something to reference', () => {
    expect(references.optionsFor('A')).toEqual([
      { node: 'C', references: [outputRef('C', 'first'), outputRef('C', 'second')] },
    ]);
  });

  it('gives what one node offers', () => {
    expect(references.referencesOf('C')).toEqual([outputRef('C', 'first'), outputRef('C', 'second')]);
    expect(references.referencesOf('B')).toEqual([]);
  });

  it('has no inline suggestions outside a reference token', () => {
    expect(references.inlineOptions('A', 'plain text', 4)).toBe(null);
  });

  it('suggests node names inside a fresh token', () => {
    expect(references.inlineOptions('A', '{', 1)).toEqual(['C']);
  });

  it('filters the node names by what was typed', () => {
    expect(references.inlineOptions('C', '{a', 2)).toEqual(['A']);
  });

  it('suggests the picked node outputs after the section dot', () => {
    expect(references.inlineOptions('A', '{C.output.', 10)).toEqual(['first', 'second']);
  });

  it('filters those outputs by what was typed', () => {
    expect(references.inlineOptions('A', '{C.output.sec', 13)).toEqual(['second']);
  });

  it('suggests nothing for a node that cannot be referenced', () => {
    expect(references.inlineOptions('A', '{B.output.', 10)).toEqual([]);
  });
});
