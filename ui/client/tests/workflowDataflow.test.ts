import { describe, it, expect } from 'vitest';
import {
  buildDataflowEdges,
  clearInputReference,
  dataflowReference,
  inputHandleId,
  insertReferenceAt,
  nodeInputHandleId,
  nodeOutputHandleId,
  outputHandleId,
  parseDataflowConnection,
  parseDataflowEdgeId,
  replaceReferenceAt,
  ReferenceTokenStage,
  setInputReference,
  tokenAtCaret,
} from '../src/components/workflow/workflowDataflow';
import { Reference } from '../src/components/workflow/references/Reference';
import { INPUT, OUTPUT } from '../src/components/workflow/references/knownKinds';
import { NodeCatalogEntry } from '../src/components/workflow/nodeCatalog';
import { NodeParameterSource, WorkflowNode } from '../src/shared/types';

// An output reference, the only kind there is today.
const outputRef = (node: string, key: string): Reference => new Reference(node, OUTPUT, key);

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
  // Kept in the legacy `parameters` form, so reading old workflows stays covered.
  A: { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: '{C.parameters.ggg}' } } },
};

describe('buildDataflowEdges', () => {
  it('links an input referencing another node output to that output', () => {
    expect(buildDataflowEdges(['C', 'A'], nodes, catalog)).toEqual([
      { id: 'df:C:out:ggg->A.bbb', source: 'C', sourceHandle: 'C:out:ggg', target: 'A', targetHandle: 'A:in:bbb' },
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

  it('ignores a reference a node makes to itself', () => {
    const n = { A: { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: '{A.output.ggg}' } } } };
    expect(buildDataflowEdges(['A'], n, catalog)).toEqual([]);
  });

  it('ignores references to a node not in the graph', () => {
    const n = { A: { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: '{Z.output.ggg}' } } } };
    expect(buildDataflowEdges(['A'], n, catalog)).toEqual([]);
  });

  it('builds one edge per reference when a value holds two', () => {
    const n = { C: nodes.C, A: { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: '{C.output.ggg} {C.output.copyDirectory}' } } } };
    expect(buildDataflowEdges(['C', 'A'], n, catalog)).toHaveLength(2);
  });

  it('builds one edge when the same reference appears twice', () => {
    const n = { C: nodes.C, A: { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: '{C.output.ggg} {C.output.ggg}' } } } };
    expect(buildDataflowEdges(['C', 'A'], n, catalog)).toHaveLength(1);
  });

  it('allows spaces inside the braces', () => {
    const n = { C: nodes.C, A: { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: '{ C.output.ggg }' } } } };
    expect(buildDataflowEdges(['C', 'A'], n, catalog)).toHaveLength(1);
  });

  it('links a reference nested inside a dict parameter (issue #1120)', () => {
    const n = { C: nodes.C, A: { type: 'general.CopyDirectory', Execution: { input_parameters: { nested: { deep: '{C.output.ggg}' } } } } };
    expect(buildDataflowEdges(['C', 'A'], n, catalog)).toEqual([
      { id: 'df:C:out:ggg->A.nested.deep', source: 'C', sourceHandle: 'C:out:ggg', target: 'A', targetHandle: 'A:in:nested.deep' },
    ]);
  });

  it('links one reference per key of a dict parameter', () => {
    const params = { Parameters: { one: '{C.output.ggg}', two: '{C.output.copyDirectory}' } };
    const n = { C: nodes.C, A: { type: 'general.CopyDirectory', Execution: { input_parameters: params } } };
    expect(buildDataflowEdges(['C', 'A'], n, catalog).map(e => e.targetHandle))
      .toEqual(['A:in:Parameters.one', 'A:in:Parameters.two']);
  });

  it('links a reference inside a list parameter, by its position', () => {
    const n = { C: nodes.C, A: { type: 'general.CopyDirectory', Execution: { input_parameters: { Command: ['echo', '{C.output.ggg}'] } } } };
    expect(buildDataflowEdges(['C', 'A'], n, catalog).map(e => e.targetHandle)).toEqual(['A:in:Command.1']);
  });

  it('ignores a parameter whose value is not a string', () => {
    const n = { C: nodes.C, A: { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: 5 } } } };
    expect(buildDataflowEdges(['C', 'A'], n, catalog)).toEqual([]);
  });

  it('has no edges for a node with no parameters at all', () => {
    expect(buildDataflowEdges(['C'], { C: nodes.C }, catalog)).toEqual([]);
  });
});

// Regression: the node-level requires handles once shared the "no id" slot with
// the dataflow handles, so a node with input/output references could no longer be
// wired with a requires edge. Every handle id must be distinct per node.
describe('handle ids let requires and dataflow coexist on one node', () => {
  it('gives a node distinct ids for its requires, output, and input handles', () => {
    const ids = [
      nodeOutputHandleId('C'),
      nodeInputHandleId('C'),
      outputHandleId('C', 'ggg'),
      inputHandleId('C', 'bbb'),
    ];
    expect(new Set(ids).size).toBe(ids.length);
  });

  it('parses a dataflow drag but not a requires drag between the same two nodes', () => {
    expect(parseDataflowConnection(nodeOutputHandleId('C'), nodeInputHandleId('A'))).toBeNull();
    expect(parseDataflowConnection(outputHandleId('C', 'ggg'), inputHandleId('A', 'bbb')))
      .toEqual({ kind: OUTPUT, outputName: 'ggg', paramPath: 'bbb' });
  });
});

describe('parseDataflowConnection', () => {
  it('parses an output→input connection into its output and param names', () => {
    expect(parseDataflowConnection('C:out:ggg', 'A:in:bbb')).toEqual({ kind: OUTPUT, outputName: 'ggg', paramPath: 'bbb' });
  });

  it('returns null when either handle is not a dataflow handle', () => {
    expect(parseDataflowConnection('C:out:ggg', null)).toBeNull();
    expect(parseDataflowConnection(null, 'A:in:bbb')).toBeNull();
    expect(parseDataflowConnection('A:in:bbb', 'C:out:ggg')).toBeNull();
  });

  it('ignores requires handles (no trailing name after the marker)', () => {
    expect(parseDataflowConnection('C:req-out', 'A:req-in')).toBeNull();
  });
});

describe('setInputReference', () => {
  it('writes into a key inside a dict parameter', () => {
    const node = { type: 'general.CopyDirectory', Execution: { input_parameters: { P: { one: '', two: 'keep' } } } };
    const updated = setInputReference(node, 'P.one', outputRef('C', 'ggg'));
    expect(updated.Execution?.input_parameters?.P).toEqual({ one: '{C.output.ggg}', two: 'keep' });
  });

  it('leaves the original node untouched', () => {
    const node = { type: 'general.CopyDirectory', Execution: { input_parameters: { P: { one: '' } } } };
    setInputReference(node, 'P.one', outputRef('C', 'ggg'));
    expect(node.Execution.input_parameters.P.one).toBe('');
  });

  it('writes {source.output.name} into the target parameter', () => {
    const updated = setInputReference({ type: 'general.CopyDirectory' }, 'bbb', outputRef('C', 'ggg'));
    expect(updated.Execution?.input_parameters?.bbb).toBe('{C.output.ggg}');
  });

  it('keeps other parameters intact', () => {
    const node = { type: 'general.CopyDirectory', Execution: { input_parameters: { aaa: '1' } } };
    const updated = setInputReference(node, 'bbb', outputRef('C', 'ggg'));
    expect(updated.Execution?.input_parameters).toEqual({ aaa: '1', bbb: '{C.output.ggg}' });
  });

  it('round-trips into a dataflow edge', () => {
    const node = setInputReference({ type: 'general.CopyDirectory' }, 'bbb', outputRef('C', 'ggg'));
    expect(buildDataflowEdges(['C', 'A'], { C: nodes.C, A: node }, catalog)).toEqual([
      { id: 'df:C:out:ggg->A.bbb', source: 'C', sourceHandle: 'C:out:ggg', target: 'A', targetHandle: 'A:in:bbb' },
    ]);
  });
});

describe('dataflowReference', () => {
  it('builds the reference token', () => {
    expect(dataflowReference('C', 'ggg')).toBe('{C.output.ggg}');
  });
});

describe('insertReferenceAt', () => {
  it('inserts the token at a caret in the middle', () => {
    expect(insertReferenceAt('ab', 1, outputRef('C', 'ggg'))).toBe('a{C.output.ggg}b');
  });

  it('inserts at the start and at the end', () => {
    expect(insertReferenceAt('ab', 0, outputRef('C', 'ggg'))).toBe('{C.output.ggg}ab');
    expect(insertReferenceAt('ab', 2, outputRef('C', 'ggg'))).toBe('ab{C.output.ggg}');
  });

  it('is just the token for an empty value', () => {
    expect(insertReferenceAt('', 0, outputRef('C', 'ggg'))).toBe('{C.output.ggg}');
  });

  it('clamps a caret out of range', () => {
    expect(insertReferenceAt('ab', -5, outputRef('C', 'ggg'))).toBe('{C.output.ggg}ab');
    expect(insertReferenceAt('ab', 99, outputRef('C', 'ggg'))).toBe('ab{C.output.ggg}');
  });
});

describe('tokenAtCaret', () => {
  it('returns null when the caret is not inside a {…} token', () => {
    expect(tokenAtCaret('hello', 3)).toBeNull();
    expect(tokenAtCaret('{C.output.ggg} tail', 16)).toBeNull();
  });

  it('reads the node stage before any dot', () => {
    expect(tokenAtCaret('{Cca', 4)).toEqual({
      stage: ReferenceTokenStage.Node, nodePart: '', sectionPart: '', kind: null, seed: 'Cca', start: 0, end: 4,
    });
  });

  it('reads the node stage right after the opening brace', () => {
    expect(tokenAtCaret('x {', 3)).toEqual({
      stage: ReferenceTokenStage.Node, nodePart: '', sectionPart: '', kind: null, seed: '', start: 2, end: 3,
    });
  });

  it('reads the section stage once a dot is typed', () => {
    expect(tokenAtCaret('{C.', 3)).toEqual({
      stage: ReferenceTokenStage.Section, nodePart: 'C', sectionPart: '', kind: null, seed: '', start: 0, end: 3,
    });
  });

  it('keeps a half-typed dotted section in the section stage', () => {
    expect(tokenAtCaret('{C.Execution.', 13)).toEqual({
      stage: ReferenceTokenStage.Section, nodePart: 'C', sectionPart: 'Execution.', kind: null, seed: 'Execution.', start: 0, end: 13,
    });
  });

  it('reads the key stage for a dotted section', () => {
    expect(tokenAtCaret('{C.Execution.input_parameters.bb', 32)).toEqual({
      stage: ReferenceTokenStage.Key, nodePart: 'C', sectionPart: 'Execution.input_parameters', kind: INPUT, seed: 'bb', start: 0, end: 32,
    });
  });

  it('filters keys by the text after the last dot', () => {
    expect(tokenAtCaret('{C.output.gg', 12)).toEqual({
      stage: ReferenceTokenStage.Key, nodePart: 'C', sectionPart: 'output', kind: OUTPUT, seed: 'gg', start: 0, end: 12,
    });
  });

  it('spans past the closing brace when the token is already closed', () => {
    const value = '{C.output.ggg}';
    expect(tokenAtCaret(value, 12)).toEqual({
      stage: ReferenceTokenStage.Key, nodePart: 'C', sectionPart: 'output', kind: OUTPUT, seed: 'gg', start: 0, end: 14,
    });
  });

  it('stops the span at the caret when the token is unclosed before another {', () => {
    expect(tokenAtCaret('{C.output.g {D', 11)).toMatchObject({ start: 0, end: 11 });
  });

  // Caret 0 sits before the brace, yet lastIndexOf searches from 0 and finds
  // it, so the node menu opens. Recorded as-is; changing it is a separate fix.
  it('treats a caret at position 0 as inside a token that starts there', () => {
    expect(tokenAtCaret('{C.output.ggg}', 0)).toMatchObject({
      stage: ReferenceTokenStage.Node, seed: '', start: 0,
    });
  });

  it('returns null at position 0 when no token starts there', () => {
    expect(tokenAtCaret('x {C.output.ggg}', 0)).toBeNull();
  });

  it('clamps a caret past the end of the value', () => {
    expect(tokenAtCaret('{C.out', 99)).toMatchObject({ nodePart: 'C', seed: 'out' });
  });
});

describe('replaceReferenceAt', () => {
  it('overwrites the token span with a full reference', () => {
    expect(replaceReferenceAt('{Cca', 0, 4, outputRef('C', 'ggg'))).toBe('{C.output.ggg}');
  });

  it('keeps text on either side of the span', () => {
    expect(replaceReferenceAt('a {C.p} b', 2, 7, outputRef('C', 'ggg'))).toBe('a {C.output.ggg} b');
  });
});

describe('parseDataflowEdgeId', () => {
  it('parses a dataflow edge id into its parts', () => {
    expect(parseDataflowEdgeId('df:C:out:ggg->A.bbb')).toEqual({ refNode: 'C', key: 'ggg', target: 'A', paramPath: 'bbb' });
  });

  it('returns null for a non-dataflow edge id', () => {
    expect(parseDataflowEdgeId('C->A')).toBeNull();
  });
});

describe('clearInputReference', () => {
  it('clears a key inside a dict parameter', () => {
    const node = { type: 'general.CopyDirectory', Execution: { input_parameters: { P: { one: '{C.output.ggg}', two: 'keep' } } } };
    const updated = clearInputReference(node, 'P.one', 'C', 'ggg');
    expect(updated.Execution?.input_parameters?.P).toEqual({ one: '', two: 'keep' });
  });

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

  it('is a no-op for a parameter that is not there', () => {
    const node = { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: 'x' } } };
    expect(clearInputReference(node, 'zzz', 'C', 'ggg').Execution?.input_parameters).toEqual({ bbb: 'x' });
  });

  it('leaves a non-string value untouched', () => {
    const node = { type: 'general.CopyDirectory', Execution: { input_parameters: { bbb: 5 } } };
    expect(clearInputReference(node, 'bbb', 'C', 'ggg').Execution?.input_parameters?.bbb).toBe(5);
  });
});
