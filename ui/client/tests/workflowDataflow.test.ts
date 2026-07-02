import { describe, it, expect } from 'vitest';
import { buildDataflowEdges } from '../src/components/workflow/workflowDataflow';
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
