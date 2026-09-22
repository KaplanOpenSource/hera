import { dataflowReference, ReferenceTokenStage, replaceReferenceAt, tokenAtCaret } from './workflowDataflow';

// What picking an inline suggestion does to a field's value: the new text, where
// the caret goes, and whether the reference is now complete (the menu closes).
export interface InlineReferencePick {
  value: string;
  caret: number;
  completed: boolean;
}

// Applies the picked suggestion to the `{…}` token the caret sits in, or returns
// null when the caret is not in one. Picking a node writes the reference
// scaffold ({node.output.}) and leaves the token open for its output; picking an
// output completes the {node.output.key} token.
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
    const scaffold = `{${option}.output.`;
    return {
      value: value.slice(0, token.start) + scaffold + value.slice(token.end),
      caret: token.start + scaffold.length,
      completed: false,
    };
  }
  return {
    value: replaceReferenceAt(value, token.start, token.end, token.nodePart, option),
    caret: token.start + dataflowReference(token.nodePart, option).length,
    completed: true,
  };
};
