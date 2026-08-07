import { describe, it, expect } from 'vitest';
import { classifyLog, LogLineKind } from '../src/components/workflow/log/classifyLog';

describe('classifyLog', () => {
  it('classifies lines by their level prefix', () => {
    const lines = classifyLog('DEBUG: checking\nINFO: running\nWARNING: hmm\nERROR: boom');
    expect(lines.map((l) => { return l.kind; })).toEqual([
      LogLineKind.Debug,
      LogLineKind.Info,
      LogLineKind.Warning,
      LogLineKind.Error,
    ]);
  });

  it('treats prefixless lines as task output', () => {
    const lines = classifyLog('hello from hera');
    expect(lines[0].kind).toBe(LogLineKind.Output);
  });

  it('indexes each line from zero', () => {
    const lines = classifyLog('a\nb\nc');
    expect(lines.map((l) => { return l.index; })).toEqual([0, 1, 2]);
  });

  it('marks the block between the two summary markers as summary', () => {
    const raw = [
      'INFO: done',
      '===== Luigi Execution Summary =====',
      'Scheduled 2 tasks',
      '===== Luigi Execution Summary =====',
      'Load /app/out.json',
    ].join('\n');
    const kinds = classifyLog(raw).map((l) => { return l.kind; });
    expect(kinds).toEqual([
      LogLineKind.Info,
      LogLineKind.Summary,
      LogLineKind.Summary,
      LogLineKind.Summary,
      LogLineKind.Output,
    ]);
  });
});
