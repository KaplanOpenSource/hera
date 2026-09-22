import { describe, it, expect } from 'vitest';
import { applyInlinePick } from '../src/components/workflow/inlineReferencePick';

describe('applyInlinePick', () => {
  it('is null when the caret is not in a reference', () => {
    expect(applyInlinePick('plain text', 5, 'A')).toBe(null);
  });

  it('writes the scaffold when a node is picked, leaving it open', () => {
    expect(applyInlinePick('{', 1, 'A')).toEqual({ value: '{A.output.', caret: 10, completed: false });
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
