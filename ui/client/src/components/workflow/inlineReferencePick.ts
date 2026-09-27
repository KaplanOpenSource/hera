import { ReferenceTokenStage, replaceReferenceAt, tokenAtCaret } from './workflowDataflow';
import { Reference } from './references/Reference';
import { knownKinds } from './references/knownKinds';

// What picking an inline suggestion does to a field's value: the new text, where
// the caret goes, and whether the reference is now complete (the menu closes).
export interface InlineReferencePick {
  value: string;
  caret: number;
  completed: boolean;
}

// Writes `text` over the token and leaves it open for the next part.
const openToken = (
  value: string,
  token: { start: number, end: number },
  text: string,
): InlineReferencePick => {
  return {
    value: value.slice(0, token.start) + text + value.slice(token.end),
    caret: token.start + text.length,
    completed: false,
  };
};

// Applies the picked suggestion to the `{…}` token the caret sits in, or returns
// null when the caret is not in one. Picking a node or a kind leaves the token
// open for the next part; picking a key completes the token.
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
    return openToken(value, token, `{${option}.`);
  }
  if (token.stage === ReferenceTokenStage.Section) {
    const kind = knownKinds.all().find(candidate => candidate.label === option);
    if (kind === undefined) {
      return null;
    }
    return openToken(value, token, `{${token.nodePart}.${kind.section}.`);
  }
  if (token.kind === null) {
    return null;
  }
  const reference = new Reference(token.nodePart, token.kind, option);
  return {
    value: replaceReferenceAt(value, token.start, token.end, reference),
    caret: token.start + reference.toString().length,
    completed: true,
  };
};
