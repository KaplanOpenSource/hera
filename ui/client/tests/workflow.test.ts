import { describe, it, expect, beforeEach } from 'vitest';
import { NO_PROJECT, useProjectStore } from '../src/stores/useProjectStore';
import { MutatorsListHandler } from '../src/shared/workflowMutators/MutatorsListHandler';
import { syncParameters, workflowParameters } from '../src/shared/workflowMutators/mutators/SyncParametersMutator';
import {
  WORKFLOW_DOC_TYPE,
  getWorkflowBlock,
  getWorkflowSolver,
  isTopLevelBlock,
  isWorkflowDoc,
  normalizeRequires,
  setWorkflowSolver,
} from '../src/shared/workflow';

describe('normalizeRequires', () => {
  it('returns an empty array for undefined', () => {
    expect(normalizeRequires(undefined)).toEqual([]);
  });

  it('wraps a single name in an array', () => {
    expect(normalizeRequires('a')).toEqual(['a']);
  });

  it('passes a list through unchanged', () => {
    expect(normalizeRequires(['a', 'b'])).toEqual(['a', 'b']);
  });
});

describe('getWorkflowBlock', () => {
  it('returns undefined for undefined', () => {
    expect(getWorkflowBlock(undefined)).toBeUndefined();
  });

  it('returns a bare block (has nodes/nodeList) as-is', () => {
    const block = { nodes: { a: {} } };
    expect(getWorkflowBlock(block)).toBe(block);
    expect(getWorkflowBlock({ nodeList: [] })).toEqual({ nodeList: [] });
  });

  it('unwraps one extra { workflow } level', () => {
    const inner = { nodes: { a: {} } };
    expect(getWorkflowBlock({ workflow: inner })).toBe(inner);
  });

  it('unwraps several nested { workflow } levels', () => {
    const inner = { nodeList: ['a'] };
    expect(getWorkflowBlock({ workflow: { workflow: inner } })).toBe(inner);
  });

  it('returns undefined when no block is present', () => {
    expect(getWorkflowBlock({})).toBeUndefined();
  });
});

describe('isTopLevelBlock', () => {
  it('is false for undefined', () => {
    expect(isTopLevelBlock(undefined)).toBe(false);
  });

  it('is true for a bare block', () => {
    expect(isTopLevelBlock({ nodes: {} })).toBe(true);
    expect(isTopLevelBlock({ nodeList: [] })).toBe(true);
  });

  it('is false for a wrapped block', () => {
    expect(isTopLevelBlock({ workflow: { nodes: {} } })).toBe(false);
  });

  it('is false for an object with neither nodes nor nodeList', () => {
    expect(isTopLevelBlock({})).toBe(false);
  });
});

describe('getWorkflowSolver', () => {
  it('returns an empty string for undefined', () => {
    expect(getWorkflowSolver(undefined)).toBe('');
  });

  it('returns the solver from a bare block', () => {
    expect(getWorkflowSolver({ nodes: {}, solver: 'simpleFoam' })).toBe('simpleFoam');
  });

  it('returns an empty string when the block has no solver', () => {
    expect(getWorkflowSolver({ nodes: {} })).toBe('');
  });

  it('reads the solver through a { workflow } wrapper', () => {
    expect(getWorkflowSolver({ workflow: { nodes: {}, solver: 'pimpleFoam' } })).toBe('pimpleFoam');
  });
});

describe('setWorkflowSolver', () => {
  it('sets the solver on a bare block, preserving other fields', () => {
    const result = setWorkflowSolver({ nodes: { a: {} }, nodeList: ['a'], solver: 'old' }, 'new');
    expect(result).toEqual({ nodes: { a: {} }, nodeList: ['a'], solver: 'new' });
  });

  it('sets the solver through a wrapper, keeping the nesting', () => {
    const result = setWorkflowSolver({ workflow: { nodes: {}, solver: 'old' } }, 'new');
    expect(result).toEqual({ workflow: { nodes: {}, solver: 'new' } });
  });

  it('produces a wrapped block when there was no workflow', () => {
    expect(setWorkflowSolver(undefined, 'new')).toEqual({ workflow: { solver: 'new' } });
  });

  it('does not mutate the input block', () => {
    const input = { nodes: {}, solver: 'old' };
    setWorkflowSolver(input, 'new');
    expect(input.solver).toBe('old');
  });
});

describe('isWorkflowDoc', () => {
  it('is false for undefined', () => {
    expect(isWorkflowDoc(undefined)).toBe(false);
  });

  it('is true when the document type matches', () => {
    expect(isWorkflowDoc({ type: WORKFLOW_DOC_TYPE })).toBe(true);
  });

  it('is true when desc.workflow holds a block', () => {
    expect(isWorkflowDoc({ desc: { workflow: { nodes: {} } } })).toBe(true);
  });

  it('is false when desc.workflow has no block', () => {
    expect(isWorkflowDoc({ desc: { workflow: {} } })).toBe(false);
  });

  it('is false for a non-workflow document', () => {
    expect(isWorkflowDoc({ type: 'other' })).toBe(false);
    expect(isWorkflowDoc({ desc: {} })).toBe(false);
  });
});

describe('workflowParameters', () => {
  const nodeWith = (parameters: { [key: string]: any }) => {
    return { Execution: { input_parameters: parameters } };
  };

  it('collects each node input_parameters keyed by node name', () => {
    const block = {
      nodeList: ['A', 'B'],
      nodes: { A: nodeWith({ ProjectName: 'P' }), B: nodeWith({ Command: 'ls' }) },
    };
    expect(workflowParameters(block)).toEqual({
      A: { ProjectName: 'P' },
      B: { Command: 'ls' },
    });
  });

  it('follows nodeList order and skips missing nodes', () => {
    const block = { nodeList: ['A', 'ghost'], nodes: { A: nodeWith({ x: 1 }) } };
    expect(workflowParameters(block)).toEqual({ A: { x: 1 } });
  });

  it('falls back to the nodes map when nodeList is absent', () => {
    const block = { nodes: { A: nodeWith({ x: 1 }) } };
    expect(workflowParameters(block)).toEqual({ A: { x: 1 } });
  });

  it('is empty for a workflow with no nodes', () => {
    expect(workflowParameters({ nodeList: [], nodes: {} })).toEqual({});
  });
});

// No project selected, so the project-name fill is a no-op and these tests see
// the parameters sync alone. The fill has its own tests in fillProjectName.
describe('MutatorsListHandler.normalize', () => {
  beforeEach(() => {
    useProjectStore.getState().selectProject(NO_PROJECT);
  });

  const descWith = (nodes: { [name: string]: any }, nodeList: string[]) => {
    return { workflow: { workflow: { nodeList, nodes } } };
  };

  it('syncs parameters (nested block)', () => {
    const desc = descWith(
      { A: { Execution: { input_parameters: { ProjectName: 'P' } } } },
      ['A'],
    );
    const result = MutatorsListHandler.normalize(desc);
    expect(result.parameters).toEqual({ A: { ProjectName: 'P' } });
  });

  it('preserves a top-level (unnested) block', () => {
    const desc = { workflow: { nodeList: ['A'], nodes: { A: { Execution: { input_parameters: { x: 1 } } } } } };
    const result = MutatorsListHandler.normalize(desc);
    expect(result.workflow?.nodeList).toEqual(['A']);
    expect(result.parameters).toEqual({ A: { x: 1 } });
  });

  it('leaves node parameters untouched when no project is selected', () => {
    const desc = descWith(
      { A: { Execution: { input_parameters: { ProjectName: '' } } } },
      ['A'],
    );
    const result = MutatorsListHandler.normalize(desc);
    expect(result.workflow?.workflow.nodes.A.Execution.input_parameters.ProjectName).toBe('');
  });

  it('returns a non-workflow desc unchanged', () => {
    const desc = { toolkit: 'X' } as any;
    expect(MutatorsListHandler.normalize(desc)).toBe(desc);
  });
});

describe('syncParameters on its own', () => {
  beforeEach(() => {
    useProjectStore.getState().selectProject(NO_PROJECT);
  });

  const desc = () => ({ workflow: { workflow: { nodeList: ['A'], nodes: { A: { Execution: { input_parameters: { ProjectName: '' } } } } } } });

  it('syncParameters only syncs parameters', () => {
    const result = syncParameters(desc());
    expect(result.workflow?.workflow.nodes.A.Execution.input_parameters.ProjectName).toBe('');
    expect(result.parameters).toEqual({ A: { ProjectName: '' } });
  });
});

// Coverage for normalize as a whole, with a project selected: both phases run,
// in order, and the parameters index is built from the FILLED values.
describe('MutatorsListHandler.normalize with a project selected', () => {
  const PROJECT = 'MY_PROJECT';

  beforeEach(() => {
    useProjectStore.getState().selectProject(PROJECT);
  });

  it('fills an empty project name and indexes the filled value', () => {
    const desc = {
      workflow: { workflow: { nodeList: ['A'], nodes: { A: { Execution: { input_parameters: { ProjectName: '' } } } } } },
    };
    const result = MutatorsListHandler.normalize(desc);
    expect(result.workflow?.workflow.nodes.A.Execution.input_parameters.ProjectName).toBe(PROJECT);
    // Proves the fill runs before the sync: the index carries the filled value.
    expect(result.parameters).toEqual({ A: { ProjectName: PROJECT } });
  });

  it('keeps a project name the user set to another project', () => {
    const desc = {
      workflow: { workflow: { nodeList: ['A'], nodes: { A: { Execution: { input_parameters: { ProjectName: 'OTHER' } } } } } },
    };
    const result = MutatorsListHandler.normalize(desc);
    expect(result.workflow?.workflow.nodes.A.Execution.input_parameters.ProjectName).toBe('OTHER');
    expect(result.parameters).toEqual({ A: { ProjectName: 'OTHER' } });
  });

  it('fills a project-name field that sits outside the workflow block', () => {
    const desc = { projectname: '', workflow: { nodeList: [], nodes: {} } };
    const result: any = MutatorsListHandler.normalize(desc);
    expect(result.projectname).toBe(PROJECT);
  });

  it('fills a project-name field nested in a list', () => {
    const desc = { steps: [{ ProjectName: '' }], workflow: { nodeList: [], nodes: {} } };
    const result: any = MutatorsListHandler.normalize(desc);
    expect(result.steps).toEqual([{ ProjectName: PROJECT }]);
  });

  it('preserves a top-level (unnested) block and other desc fields', () => {
    const desc = {
      workflowName: 'W',
      workflow: { nodeList: ['A'], nodes: { A: { Execution: { input_parameters: { x: 1 } } } } },
    };
    const result: any = MutatorsListHandler.normalize(desc);
    expect(result.workflowName).toBe('W');
    expect(result.workflow.nodeList).toEqual(['A']);
    expect(result.parameters).toEqual({ A: { x: 1 } });
  });

  it('handles a block with no nodes', () => {
    const desc = { workflow: { nodeList: [], nodes: {} } };
    const result: any = MutatorsListHandler.normalize(desc);
    expect(result.parameters).toEqual({});
  });

  it('leaves a node that has no Execution out of the index', () => {
    const desc = { workflow: { nodeList: ['A'], nodes: { A: { type: 'general.RunOsCommand' } } } };
    const result: any = MutatorsListHandler.normalize(desc);
    // toStrictEqual, so an { A: undefined } entry is a failure, not a match.
    expect(result.parameters).toStrictEqual({});
    expect(result.workflow.nodes.A).toEqual({ type: 'general.RunOsCommand' });
  });

  it('indexes only the nodes in nodeList, in that order', () => {
    const desc = {
      workflow: {
        nodeList: ['B', 'A'],
        nodes: {
          A: { Execution: { input_parameters: { x: 1 } } },
          B: { Execution: { input_parameters: { y: 2 } } },
          C: { Execution: { input_parameters: { z: 3 } } },
        },
      },
    };
    const result: any = MutatorsListHandler.normalize(desc);
    // C is not in nodeList, so it is not indexed; B comes before A.
    expect(Object.keys(result.parameters)).toEqual(['B', 'A']);
  });

  it('returns a non-workflow desc unchanged', () => {
    const desc = { toolkit: 'X' } as any;
    const result: any = MutatorsListHandler.normalize(desc);
    expect(result).toEqual({ toolkit: 'X' });
    expect(result.parameters).toBeUndefined();
  });

  it('fills a project name that is null or unset', () => {
    const desc = {
      workflow: { nodeList: ['A', 'B'], nodes: {
        A: { Execution: { input_parameters: { ProjectName: null } } },
        B: { Execution: { input_parameters: { ProjectName: undefined } } },
      } },
    };
    const result: any = MutatorsListHandler.normalize(desc);
    expect(result.workflow.nodes.A.Execution.input_parameters.ProjectName).toBe(PROJECT);
    expect(result.workflow.nodes.B.Execution.input_parameters.ProjectName).toBe(PROJECT);
  });

  it('matches the field name in any casing, and only that name', () => {
    const desc = { PROJECTNAME: '', projectNames: '', theProjectName: '', workflow: { nodeList: [], nodes: {} } };
    const result: any = MutatorsListHandler.normalize(desc);
    expect(result.PROJECTNAME).toBe(PROJECT);
    // Near misses stay empty: the name must match end to end.
    expect(result.projectNames).toBe('');
    expect(result.theProjectName).toBe('');
  });

  it('leaves a null value under another field alone', () => {
    const desc = { note: null, workflow: { nodeList: [], nodes: {} } };
    const result: any = MutatorsListHandler.normalize(desc);
    expect(result.note).toBeNull();
  });

  it('handles a block that has nodeList but no nodes map', () => {
    const desc = { workflow: { nodeList: ['A'] } };
    const result: any = MutatorsListHandler.normalize(desc);
    expect(result.parameters).toStrictEqual({});
  });

  it('falls back to the nodes map when the block has no nodeList', () => {
    const desc = { workflow: { nodes: { A: { Execution: { input_parameters: { x: 1 } } } } } };
    const result: any = MutatorsListHandler.normalize(desc);
    expect(result.parameters).toEqual({ A: { x: 1 } });
  });

  it('does not mutate the desc it was given', () => {
    const desc = {
      workflow: { workflow: { nodeList: ['A'], nodes: { A: { Execution: { input_parameters: { ProjectName: '' } } } } } },
    };
    MutatorsListHandler.normalize(desc);
    expect(desc.workflow.workflow.nodes.A.Execution.input_parameters.ProjectName).toBe('');
    expect((desc as any).parameters).toBeUndefined();
  });
});

// The remaining guard in the fill phase: an empty project name in the store.
describe('MutatorsListHandler.normalize with an empty project name', () => {
  beforeEach(() => {
    useProjectStore.getState().selectProject('');
  });

  it('does not fill, but still syncs parameters', () => {
    const desc = { workflow: { nodeList: ['A'], nodes: { A: { Execution: { input_parameters: { ProjectName: '' } } } } } };
    const result: any = MutatorsListHandler.normalize(desc);
    expect(result.workflow.nodes.A.Execution.input_parameters.ProjectName).toBe('');
    expect(result.parameters).toEqual({ A: { ProjectName: '' } });
  });
});
