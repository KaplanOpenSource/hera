import { describe, it, expect } from 'vitest';
import { applyInlinePick } from '../src/components/workflow/inlineReferencePick';

describe('applyInlinePick', () => {
  it('is null when the caret is not in a reference', () => {
    expect(applyInlinePick('plain text', 5, 'A')).toBe(null);
  });

  it('writes the node when a node is picked, leaving it open', () => {
    expect(applyInlinePick('{', 1, 'A')).toEqual({ value: '{A.', caret: 3, completed: false });
  });

  it('writes the section when a kind is picked, leaving it open', () => {
    expect(applyInlinePick('{A.', 3, 'Output')).toEqual({ value: '{A.output.', caret: 10, completed: false });
    expect(applyInlinePick('{A.', 3, 'Input')).toEqual({
      value: '{A.Execution.input_parameters.',
      caret: 30,
      completed: false,
    });
  });

  it('completes an input reference when its key is picked', () => {
    expect(applyInlinePick('{A.Execution.input_parameters.', 30, 'cmd')).toEqual({
      value: '{A.Execution.input_parameters.cmd}',
      caret: 34,
      completed: true,
    });
  });

  it('completes the reference when an output is picked', () => {
    expect(applyInlinePick('{A.output.', 10, 'result')).toEqual({
      value: '{A.output.result}',
      caret: 17,
      completed: true,
    });
  });

  it('keeps the text around the reference', () => {
    expect(applyInlinePick('cmd {A.output. --flag', 14, 'result')).toEqual({
      value: 'cmd {A.output.result} --flag',
      caret: 21,
      completed: true,
    });
  });

  it('replaces a half-typed output', () => {
    expect(applyInlinePick('{A.output.res', 13, 'result').value).toBe('{A.output.result}');
  });
});
