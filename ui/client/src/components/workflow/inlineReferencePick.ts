import { ReferenceTokenStage, replaceReferenceAt, tokenAtCaret } from './workflowDataflow';
import { Reference } from './references/Reference';
import { OUTPUT } from './references/knownKinds';

// What picking an inline suggestion does to a field's value: the new text, where
// the caret goes, and whether the reference is now complete (the menu closes).
export interface InlineReferencePick {
  value: string;
  caret: number;
  completed: boolean;
}

// Applies the picked suggestion to the `{…}` token the caret sits in, or returns
// null when the caret is not in one. Picking a node writes the reference scaffold
// and leaves the token open for its key; picking a key completes the token.
export const applyInlinePick = (
  value: string,
  caret: number,
  option: string,
): InlineReferencePick | null => {
  const token = tokenAtCaret(value, caret);
  if (token === null) {
    return null;
  }
  if (token.stage === ReferenceTokenStage.Node) {
    const scaffold = `{${option}.${OUTPUT.section}.`;
    return {
      value: value.slice(0, token.start) + scaffold + value.slice(token.end),
      caret: token.start + scaffold.length,
      completed: false,
    };
  }
  const reference = new Reference(token.nodePart, OUTPUT, option);
  return {
    value: replaceReferenceAt(value, token.start, token.end, reference),
    caret: token.start + reference.toString().length,
    completed: true,
  };
};
