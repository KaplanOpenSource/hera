import { describe, it, expect } from 'vitest';
import {
  WORKFLOW_DOC_TYPE,
  fillNewProjectNameParams,
  fillProjectName,
  isProjectNameKey,
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

describe('isProjectNameKey', () => {
  it('matches ProjectName in any casing', () => {
    expect(isProjectNameKey('ProjectName')).toBe(true);
    expect(isProjectNameKey('projectName')).toBe(true);
    expect(isProjectNameKey('projectname')).toBe(true);
    expect(isProjectNameKey('PROJECTNAME')).toBe(true);
  });

  it('does not match other keys', () => {
    expect(isProjectNameKey('project')).toBe(false);
    expect(isProjectNameKey('projectNames')).toBe(false);
    expect(isProjectNameKey('SimulationName')).toBe(false);
  });
});

describe('fillProjectName', () => {
  const nodeWith = (parameters: { [key: string]: any }) => {
    return { Execution: { input_parameters: parameters } };
  };

  it('fills a ProjectName param that was empty', () => {
    const block = { nodes: { A: nodeWith({ ProjectName: '' }) } };
    const result = fillProjectName(block, 'MY_PROJECT');
    expect(result.nodes?.A.Execution?.input_parameters?.ProjectName).toBe('MY_PROJECT');
  });

  it('overwrites an existing ProjectName value', () => {
    const block = { nodes: { A: nodeWith({ ProjectName: 'OTHER' }) } };
    const result = fillProjectName(block, 'MY_PROJECT');
    expect(result.nodes?.A.Execution?.input_parameters?.ProjectName).toBe('MY_PROJECT');
  });

  it('matches the key in any casing', () => {
    const block = { nodes: { A: nodeWith({ projectName: '' }) } };
    const result = fillProjectName(block, 'MY_PROJECT');
    expect(result.nodes?.A.Execution?.input_parameters?.projectName).toBe('MY_PROJECT');
  });

  it('preserves other parameters', () => {
    const block = { nodes: { A: nodeWith({ ProjectName: '', Template: 'T' }) } };
    const result = fillProjectName(block, 'MY_PROJECT');
    expect(result.nodes?.A.Execution?.input_parameters).toEqual({ ProjectName: 'MY_PROJECT', Template: 'T' });
  });

  it('leaves a node with no project param untouched', () => {
    const block = { nodes: { A: nodeWith({ Command: 'ls' }) } };
    const result = fillProjectName(block, 'MY_PROJECT');
    expect(result.nodes?.A.Execution?.input_parameters).toEqual({ Command: 'ls' });
  });

  it('leaves a node without input_parameters untouched', () => {
    const block = { nodes: { A: { type: 'general.RunOsCommand' } } };
    const result = fillProjectName(block, 'MY_PROJECT');
    expect(result.nodes?.A).toEqual({ type: 'general.RunOsCommand' });
  });

  it('fills every matching node', () => {
    const block = { nodes: { A: nodeWith({ ProjectName: '' }), B: nodeWith({ projectName: 'OLD' }) } };
    const result = fillProjectName(block, 'MY_PROJECT');
    expect(result.nodes?.A.Execution?.input_parameters?.ProjectName).toBe('MY_PROJECT');
    expect(result.nodes?.B.Execution?.input_parameters?.projectName).toBe('MY_PROJECT');
  });

  it('does not mutate the input block', () => {
    const params = { ProjectName: '' };
    const block = { nodes: { A: nodeWith(params) } };
    fillProjectName(block, 'MY_PROJECT');
    expect(params).toEqual({ ProjectName: '' });
  });
});

describe('fillNewProjectNameParams', () => {
  it('fills a project param that was just created', () => {
    const result = fillNewProjectNameParams({ ProjectName: '' }, {}, 'MY_PROJECT');
    expect(result).toEqual({ ProjectName: 'MY_PROJECT' });
  });

  it('fills a project param renamed into existence', () => {
    const result = fillNewProjectNameParams({ ProjectName: '' }, { foo: 'x' }, 'MY_PROJECT');
    expect(result).toEqual({ ProjectName: 'MY_PROJECT' });
  });

  it('leaves an existing project param alone so the user can override it', () => {
    const result = fillNewProjectNameParams({ ProjectName: 'TYPED' }, { ProjectName: '' }, 'MY_PROJECT');
    expect(result).toEqual({ ProjectName: 'TYPED' });
  });

  it('matches the key in any casing', () => {
    const result = fillNewProjectNameParams({ projectname: '' }, {}, 'MY_PROJECT');
    expect(result).toEqual({ projectname: 'MY_PROJECT' });
  });

  it('returns the same object when nothing changed', () => {
    const next = { Command: 'ls' };
    expect(fillNewProjectNameParams(next, {}, 'MY_PROJECT')).toBe(next);
  });
});
