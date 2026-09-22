import { describe, it, expect } from 'vitest';
import { nodeNameFromTask, taskBelongsToNode } from '../src/components/workflow/taskNodeName';

const nodeNames = ['ListFiles', 'finalnode_xx', 'Run', 'Run_extra'];

describe('nodeNameFromTask', () => {
  it('maps a generated task back to its node', () => {
    expect(nodeNameFromTask('ListFiles_0', nodeNames)).toBe('ListFiles');
  });

  it('keeps a task name that is already a node name', () => {
    expect(nodeNameFromTask('ListFiles', nodeNames)).toBe('ListFiles');
  });

  it('handles a node name that itself contains an underscore', () => {
    expect(nodeNameFromTask('finalnode_xx_0', nodeNames)).toBe('finalnode_xx');
  });

  it('prefers the longest matching node', () => {
    expect(nodeNameFromTask('Run_extra_0', nodeNames)).toBe('Run_extra');
  });

  it('finds nothing for a task that is not a node', () => {
    expect(nodeNameFromTask('__between__', nodeNames)).toBeUndefined();
    expect(nodeNameFromTask('Missing_0', nodeNames)).toBeUndefined();
  });

  it('does not match a suffix that is not an index', () => {
    expect(nodeNameFromTask('ListFiles_final', nodeNames)).toBeUndefined();
  });
});

describe('taskBelongsToNode', () => {
  it('matches a generated task to its node', () => {
    expect(taskBelongsToNode('ListFiles_0', 'ListFiles')).toBe(true);
  });

  it('matches a task named exactly like the node', () => {
    expect(taskBelongsToNode('ListFiles', 'ListFiles')).toBe(true);
  });

  it('does not match another node that shares a prefix', () => {
    expect(taskBelongsToNode('Run_extra_0', 'Run')).toBe(false);
  });

  it('does not match filler output', () => {
    expect(taskBelongsToNode('__between__', 'ListFiles')).toBe(false);
  });
});
