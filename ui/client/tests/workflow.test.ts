import { describe, it, expect, beforeEach } from 'vitest';
import { NO_PROJECT, useProjectStore } from '../src/stores/useProjectStore';
import { MutatorsListHandler } from '../src/shared/workflowMutators/MutatorsListHandler';
import { WorkflowMutatorBase } from '../src/shared/workflowMutators/WorkflowMutatorBase';
import { SyncParametersMutator, workflowParameters } from '../src/shared/workflowMutators/mutators/SyncParametersMutator';
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

describe('workflow mutator phases', () => {
  beforeEach(() => {
    useProjectStore.getState().selectProject(NO_PROJECT);
  });

  const desc = () => ({ workflow: { workflow: { nodeList: ['A'], nodes: { A: { Execution: { input_parameters: { ProjectName: '' } } } } } } });

  it('exposes each phase as a WorkflowMutator', () => {
    expect(MutatorsListHandler.mutators.every(m => m instanceof WorkflowMutatorBase)).toBe(true);
    expect(MutatorsListHandler.mutators.map(m => m.name)).toEqual(['fillProjectName', 'syncParameters']);
  });

  it('SyncParametersMutator only syncs parameters', () => {
    const result = new SyncParametersMutator().mutate(desc());
    expect(result.workflow?.workflow.nodes.A.Execution.input_parameters.ProjectName).toBe('');
    expect(result.parameters).toEqual({ A: { ProjectName: '' } });
  });
});
