import { describe, it, expect } from 'vitest';
import { nodeRunStatuses, NodeRunStatus } from '../src/components/workflow/nodeRunStatus';

const nodeNames = ['ListFiles', 'Run', 'Run_extra'];

const chunk = (name: string, text: string) => {
  return { name, text };
};

describe('nodeRunStatuses', () => {
  it('has no status for a node with no events', () => {
    expect(nodeRunStatuses({ chunks: [], nodeNames, runFinished: false })).toEqual({});
  });

  it('marks a started task as running', () => {
    const chunks = [chunk('ListFiles_0', '[luigi-event] START ListFiles_0\nhello\n')];
    expect(nodeRunStatuses({ chunks, nodeNames, runFinished: false })).toEqual({
      ListFiles: NodeRunStatus.Running,
    });
  });

  it('marks a finished task as success', () => {
    const chunks = [chunk('ListFiles_0', '[luigi-event] START ListFiles_0\n[luigi-event] SUCCESS ListFiles_0\n')];
    expect(nodeRunStatuses({ chunks, nodeNames, runFinished: true })).toEqual({
      ListFiles: NodeRunStatus.Success,
    });
  });

  it('marks a failed task as failure, with its message on the line', () => {
    const chunks = [chunk('Run_0', '[luigi-event] START Run_0\n[luigi-event] FAILURE Run_0: boom\n')];
    expect(nodeRunStatuses({ chunks, nodeNames, runFinished: true })).toEqual({
      Run: NodeRunStatus.Failure,
    });
  });

  it('treats a broken task as failure', () => {
    const chunks = [chunk('Run_0', '[luigi-event] BROKEN Run_0: bad dependency\n')];
    expect(nodeRunStatuses({ chunks, nodeNames, runFinished: true })).toEqual({
      Run: NodeRunStatus.Failure,
    });
  });

  it('maps the longest matching node name', () => {
    const chunks = [chunk('Run_extra_0', '[luigi-event] SUCCESS Run_extra_0\n')];
    expect(nodeRunStatuses({ chunks, nodeNames, runFinished: true })).toEqual({
      Run_extra: NodeRunStatus.Success,
    });
  });

  it('keeps each node of several tasks apart', () => {
    const chunks = [
      chunk('ListFiles_0', '[luigi-event] START ListFiles_0\n[luigi-event] SUCCESS ListFiles_0\n'),
      chunk('Run_0', '[luigi-event] START Run_0\n'),
    ];
    expect(nodeRunStatuses({ chunks, nodeNames, runFinished: false })).toEqual({
      ListFiles: NodeRunStatus.Success,
      Run: NodeRunStatus.Running,
    });
  });

  it('counts a task still running when the run ended as failed', () => {
    const chunks = [chunk('Run_0', '[luigi-event] START Run_0\n')];
    expect(nodeRunStatuses({ chunks, nodeNames, runFinished: true })).toEqual({
      Run: NodeRunStatus.Failure,
    });
  });

  it('ignores lines that are not luigi events', () => {
    const chunks = [chunk('Run_0', 'START Run_0\nINFO: SUCCESS Run_0\n')];
    expect(nodeRunStatuses({ chunks, nodeNames, runFinished: false })).toEqual({});
  });
});
