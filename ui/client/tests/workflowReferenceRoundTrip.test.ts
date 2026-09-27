import { describe, it, expect } from 'vitest';
import {
  buildDataflowEdges,
  clearInputReference,
  dataflowReference,
  inputHandleId,
  insertReferenceAt,
  outputHandleId,
  parseDataflowConnection,
  parseDataflowEdgeId,
  setInputReference,
  tokenAtCaret,
} from '../src/components/workflow/workflowDataflow';
import { applyInlinePick } from '../src/components/workflow/inlineReferencePick';
import { WorkflowReferences } from '../src/components/workflow/WorkflowReferences';
import { Reference } from '../src/components/workflow/references/Reference';
import { OUTPUT } from '../src/components/workflow/references/knownKinds';
import { NodeCatalogEntry } from '../src/components/workflow/nodeCatalog';
import { NodeParameterSource, WorkflowNode } from '../src/shared/types';

// A reference is written in one place, read back in five others: the value
// parser, the handle ids, the edge id, the clear helper, and the menus. These
// tests pin those five against the writer without naming any of the formats,
// so the formats stay free to change.

const catalog: NodeCatalogEntry[] = [
  {
    type: 'maker',
    parameters: [],
    outputs: [
      { name: 'result', source: NodeParameterSource.Python },
      { name: 'other', source: NodeParameterSource.Python },
    ],
  },
  { type: 'plain', parameters: [] },
];

const NAMES = ['A', 'B'];

// An output reference, the only kind there is today.
const outputRef = (node: string, key: string): Reference => new Reference(node, OUTPUT, key);

// Node B with one parameter, and node A producing the outputs above.
const workflow = (cmd: any): { [name: string]: WorkflowNode } => {
  return {
    A: { type: 'maker' },
    B: { type: 'plain', Execution: { input_parameters: { cmd } } },
  };
};

// The single edge the workflow should produce, failing loudly when it does not.
const onlyEdge = (nodes: { [name: string]: WorkflowNode }) => {
  const edges = buildDataflowEdges(NAMES, nodes, catalog);
  expect(edges).toHaveLength(1);
  return edges[0];
};

describe('a written reference is found by the value parser', () => {
  it('links the node that wrote it to the node it names', () => {
    const edge = onlyEdge(workflow(dataflowReference('A', 'result')));
    expect(edge.source).toBe('A');
    expect(edge.target).toBe('B');
  });

  it('is found with text around it', () => {
    const edge = onlyEdge(workflow(`run ${dataflowReference('A', 'result')} now`));
    expect(edge.source).toBe('A');
  });

  it('is found through setInputReference too', () => {
    const nodes = workflow('');
    nodes.B = setInputReference(nodes.B, 'cmd', outputRef('A', 'result'));
    expect(onlyEdge(nodes).source).toBe('A');
  });
});

describe('an edge id survives a round trip', () => {
  it('parses back into the reference it was built from', () => {
    const edge = onlyEdge(workflow(dataflowReference('A', 'result')));
    expect(parseDataflowEdgeId(edge.id)).toMatchObject({
      refNode: 'A',
      key: 'result',
      target: 'B',
      param: 'cmd',
    });
  });

  it('is not mistaken for a requires edge id', () => {
    expect(parseDataflowEdgeId('A->B')).toBeNull();
  });
});

describe('an edge hangs off the same dots a drag uses', () => {
  it('reads its own handles back as a dataflow connection', () => {
    const edge = onlyEdge(workflow(dataflowReference('A', 'result')));
    expect(parseDataflowConnection(edge.sourceHandle, edge.targetHandle)).toEqual({
      kind: OUTPUT,
      outputName: 'result',
      param: 'cmd',
    });
  });

  it('puts the edge on the handles the row components render', () => {
    const edge = onlyEdge(workflow(dataflowReference('A', 'result')));
    expect(edge.sourceHandle).toBe(outputHandleId('A', 'result'));
    expect(edge.targetHandle).toBe(inputHandleId('B', 'cmd'));
  });
});

describe('dragging a line and drawing it agree', () => {
  it('turns a drag between two handles into the very same edge', () => {
    const drag = parseDataflowConnection(outputHandleId('A', 'result'), inputHandleId('B', 'cmd'));
    const nodes = workflow('');
    nodes.B = setInputReference(nodes.B, drag!.param, new Reference('A', drag!.kind, drag!.outputName));
    const edge = onlyEdge(nodes);
    expect(edge.sourceHandle).toBe(outputHandleId('A', 'result'));
    expect(edge.targetHandle).toBe(inputHandleId('B', 'cmd'));
  });
});

describe('deleting a line undoes what writing it did', () => {
  it('clears the value and the edge with it', () => {
    const nodes = workflow('');
    nodes.B = setInputReference(nodes.B, 'cmd', outputRef('A', 'result'));
    const parsed = parseDataflowEdgeId(onlyEdge(nodes).id)!;
    nodes.B = clearInputReference(nodes.B, parsed.param, parsed.refNode, parsed.key);
    expect(nodes.B.Execution?.input_parameters?.cmd).toBe('');
    expect(buildDataflowEdges(NAMES, nodes, catalog)).toEqual([]);
  });

  it('leaves the surrounding text behind', () => {
    const nodes = workflow(`run ${dataflowReference('A', 'result')} now`);
    const parsed = parseDataflowEdgeId(onlyEdge(nodes).id)!;
    nodes.B = clearInputReference(nodes.B, parsed.param, parsed.refNode, parsed.key);
    expect(nodes.B.Execution?.input_parameters?.cmd).toBe('run  now');
  });

  it('leaves a second reference in the same value alone', () => {
    const value = `${dataflowReference('A', 'result')} ${dataflowReference('A', 'other')}`;
    const nodes = workflow(value);
    nodes.B = clearInputReference(nodes.B, 'cmd', 'A', 'result');
    expect(buildDataflowEdges(NAMES, nodes, catalog)).toHaveLength(1);
  });
});

describe('the typed token and the written token line up', () => {
  it('sees the caret as inside a token just inserted', () => {
    const value = insertReferenceAt('', 0, outputRef('A', 'result'));
    const token = tokenAtCaret(value, value.length - 1);
    expect(token).not.toBeNull();
    expect(token!.nodePart).toBe('A');
  });

  it('spans the whole token it just wrote', () => {
    const value = insertReferenceAt('xx', 2, outputRef('A', 'result'));
    const token = tokenAtCaret(value, value.length - 1)!;
    expect(value.slice(token.start, token.end)).toBe(dataflowReference('A', 'result'));
  });
});

describe('picking through the inline menu writes a resolvable reference', () => {
  it('builds an edge from a node pick, a kind pick then an output pick', () => {
    const afterNode = applyInlinePick('{', 1, 'A')!;
    expect(afterNode.completed).toBe(false);
    const afterKind = applyInlinePick(afterNode.value, afterNode.caret, 'Output')!;
    expect(afterKind.completed).toBe(false);
    const afterOutput = applyInlinePick(afterKind.value, afterKind.caret, 'result')!;
    expect(afterOutput.completed).toBe(true);
    expect(onlyEdge(workflow(afterOutput.value)).source).toBe('A');
  });

  it('leaves the caret past the finished token', () => {
    const afterNode = applyInlinePick('{', 1, 'A')!;
    const afterKind = applyInlinePick(afterNode.value, afterNode.caret, 'Output')!;
    const afterOutput = applyInlinePick(afterKind.value, afterKind.caret, 'result')!;
    expect(afterOutput.caret).toBe(afterOutput.value.length);
  });
});

describe('everything the menus offer can actually be drawn', () => {
  const references = new WorkflowReferences(NAMES, workflow(''), catalog);

  it('turns every offered reference into a real edge', () => {
    const offered = references.optionsFor('B');
    expect(offered.length).toBeGreaterThan(0);
    offered.forEach(option => {
      option.references.forEach(offered => {
        expect(buildDataflowEdges(NAMES, workflow(offered.toString()), catalog)).toHaveLength(1);
      });
    });
  });

  it('never offers the node being edited', () => {
    expect(references.optionsFor('A').map(option => option.node)).not.toContain('A');
  });
});
