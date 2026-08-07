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

  it('treats CRITICAL as error and WARN as warning', () => {
    const lines = classifyLog('CRITICAL: fatal\nWARN: heads up');
    expect(lines.map((l) => { return l.kind; })).toEqual([
      LogLineKind.Error,
      LogLineKind.Warning,
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
      LogLineKind.Technical,
    ]);
  });

  it('classifies everything after the summary block as post-summary', () => {
    const raw = [
      '===== Luigi Execution Summary =====',
      '===== Luigi Execution Summary =====',
      'Load /app/out.json',
      'INFO: this is after the summary too',
    ].join('\n');
    const kinds = classifyLog(raw).map((l) => { return l.kind; });
    expect(kinds).toEqual([
      LogLineKind.Summary,
      LogLineKind.Summary,
      LogLineKind.Technical,
      LogLineKind.Technical,
    ]);
  });
});
