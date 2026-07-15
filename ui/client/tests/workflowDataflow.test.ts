import { describe, it, expect } from 'vitest';
import {
  buildDataflowEdges,
  clearInputReference,
  dataflowReference,
  insertReferenceAt,
  parseDataflowConnection,
  parseDataflowEdgeId,
  replaceReferenceAt,
  ReferenceTokenStage,
  setInputReference,
  tokenAtCaret,
} from '../src/components/workflow/workflowDataflow';
import { NodeCatalogEntry } from '../src/components/workflow/nodeCatalog';
import { NodeParameterSource, WorkflowNode } from '../src/shared/types';

const catalog: NodeCatalogEntry[] = [{
  type: 'general.CopyDirectory',
  parameters: [],
  outputs: [
    { name: 'ggg', source: NodeParameterSource.Python },
    { name: 'copyDirectory', source: NodeParameterSource.Python },
  ],
}];

const nodes: { [name: string]: WorkflowNode } = {
  C: { type: 'general.CopyDirectory' },
  A: { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: '{C.parameters.ggg}' } } },
};

describe('buildDataflowEdges', () => {
  it('links an input referencing another node output to that output', () => {
    expect(buildDataflowEdges(['C', 'A'], nodes, catalog)).toEqual([
      { id: 'df:C.ggg->A.bbb', source: 'C', sourceHandle: 'out:ggg', target: 'A', targetHandle: 'in:bbb' },
    ]);
  });

  it('matches the reference embedded in a larger string', () => {
    const n = { A: { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: 'x {C.output.ggg} y' } } }, C: nodes.C };
    expect(buildDataflowEdges(['C', 'A'], n, catalog)).toHaveLength(1);
  });

  it('ignores references to a key that is not an output', () => {
    const n = { A: { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: '{C.output.nope}' } } }, C: nodes.C };
    expect(buildDataflowEdges(['C', 'A'], n, catalog)).toEqual([]);
  });

  it('ignores references to a node not in the graph', () => {
    const n = { A: { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: '{Z.output.ggg}' } } } };
    expect(buildDataflowEdges(['A'], n, catalog)).toEqual([]);
  });
});

describe('parseDataflowConnection', () => {
  it('parses an output→input connection into its output and param names', () => {
    expect(parseDataflowConnection('out:ggg', 'in:bbb')).toEqual({ outputName: 'ggg', param: 'bbb' });
  });

  it('returns null when either handle is not a dataflow handle', () => {
    expect(parseDataflowConnection('out:ggg', null)).toBeNull();
    expect(parseDataflowConnection(null, 'in:bbb')).toBeNull();
    expect(parseDataflowConnection('in:bbb', 'out:ggg')).toBeNull();
  });
});

describe('setInputReference', () => {
  it('writes {source.parameters.name} into the target parameter', () => {
    const updated = setInputReference({ type: 'general.CopyDirectory' }, 'bbb', 'C', 'ggg');
    expect(updated.Execution?.input_parameters?.bbb).toBe('{C.parameters.ggg}');
  });

  it('keeps other parameters intact', () => {
    const node = { type: 'general.CopyDirectory', Execution: { input_parameters: { aaa: '1' } } };
    const updated = setInputReference(node, 'bbb', 'C', 'ggg');
    expect(updated.Execution?.input_parameters).toEqual({ aaa: '1', bbb: '{C.parameters.ggg}' });
  });

  it('round-trips into a dataflow edge', () => {
    const node = setInputReference({ type: 'general.CopyDirectory' }, 'bbb', 'C', 'ggg');
    expect(buildDataflowEdges(['C', 'A'], { C: nodes.C, A: node }, catalog)).toEqual([
      { id: 'df:C.ggg->A.bbb', source: 'C', sourceHandle: 'out:ggg', target: 'A', targetHandle: 'in:bbb' },
    ]);
  });
});

describe('dataflowReference', () => {
  it('builds the reference token', () => {
    expect(dataflowReference('C', 'ggg')).toBe('{C.parameters.ggg}');
  });
});

describe('insertReferenceAt', () => {
  it('inserts the token at a caret in the middle', () => {
    expect(insertReferenceAt('ab', 1, 'C', 'ggg')).toBe('a{C.parameters.ggg}b');
  });

  it('inserts at the start and at the end', () => {
    expect(insertReferenceAt('ab', 0, 'C', 'ggg')).toBe('{C.parameters.ggg}ab');
    expect(insertReferenceAt('ab', 2, 'C', 'ggg')).toBe('ab{C.parameters.ggg}');
  });

  it('is just the token for an empty value', () => {
    expect(insertReferenceAt('', 0, 'C', 'ggg')).toBe('{C.parameters.ggg}');
  });

  it('clamps a caret out of range', () => {
    expect(insertReferenceAt('ab', -5, 'C', 'ggg')).toBe('{C.parameters.ggg}ab');
    expect(insertReferenceAt('ab', 99, 'C', 'ggg')).toBe('ab{C.parameters.ggg}');
  });
});

describe('tokenAtCaret', () => {
  it('returns null when the caret is not inside a {…} token', () => {
    expect(tokenAtCaret('hello', 3)).toBeNull();
    expect(tokenAtCaret('{C.parameters.ggg} tail', 20)).toBeNull();
  });

  it('reads the node stage before any dot', () => {
    expect(tokenAtCaret('{Cca', 4)).toEqual({
      stage: ReferenceTokenStage.Node, nodePart: '', seed: 'Cca', start: 0, end: 4,
    });
  });

  it('reads the node stage right after the opening brace', () => {
    expect(tokenAtCaret('x {', 3)).toEqual({
      stage: ReferenceTokenStage.Node, nodePart: '', seed: '', start: 2, end: 3,
    });
  });

  it('reads the output stage once a dot is typed', () => {
    expect(tokenAtCaret('{C.', 3)).toEqual({
      stage: ReferenceTokenStage.Output, nodePart: 'C', seed: '', start: 0, end: 3,
    });
  });

  it('filters output keys by the text after the last dot', () => {
    expect(tokenAtCaret('{C.parameters.gg', 16)).toEqual({
      stage: ReferenceTokenStage.Output, nodePart: 'C', seed: 'gg', start: 0, end: 16,
    });
  });

  it('spans past the closing brace when the token is already closed', () => {
    const value = '{C.parameters.ggg}';
    expect(tokenAtCaret(value, 16)).toEqual({
      stage: ReferenceTokenStage.Output, nodePart: 'C', seed: 'gg', start: 0, end: 18,
    });
  });

  it('stops the span at the caret when the token is unclosed before another {', () => {
    expect(tokenAtCaret('{C.parameters.g {D', 15)).toMatchObject({ start: 0, end: 15 });
  });
});

describe('replaceReferenceAt', () => {
  it('overwrites the token span with a full reference', () => {
    expect(replaceReferenceAt('{Cca', 0, 4, 'C', 'ggg')).toBe('{C.parameters.ggg}');
  });

  it('keeps text on either side of the span', () => {
    expect(replaceReferenceAt('a {C.p} b', 2, 7, 'C', 'ggg')).toBe('a {C.parameters.ggg} b');
  });
});

describe('parseDataflowEdgeId', () => {
  it('parses a dataflow edge id into its parts', () => {
    expect(parseDataflowEdgeId('df:C.ggg->A.bbb')).toEqual({ refNode: 'C', key: 'ggg', target: 'A', param: 'bbb' });
  });

  it('returns null for a non-dataflow edge id', () => {
    expect(parseDataflowEdgeId('C->A')).toBeNull();
  });
});

describe('clearInputReference', () => {
  it('clears a parameter that is exactly the reference', () => {
    const node = { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: '{C.output.ggg}' } } };
    const updated = clearInputReference(node, 'bbb', 'C', 'ggg');
    expect(updated.Execution?.input_parameters?.bbb).toBe('');
  });

  it('strips only the matching reference from a larger string', () => {
    const node = { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: 'x {C.output.ggg} y' } } };
    const updated = clearInputReference(node, 'bbb', 'C', 'ggg');
    expect(updated.Execution?.input_parameters?.bbb).toBe('x  y');
  });

  it('leaves a non-string value untouched', () => {
    const node = { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: 5 } } };
    expect(clearInputReference(node, 'bbb', 'C', 'ggg').Execution?.input_parameters?.bbb).toBe(5);
  });
});
