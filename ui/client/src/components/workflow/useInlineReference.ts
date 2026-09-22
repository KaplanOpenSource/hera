import { useEffect, useRef, useState } from 'react';
import { WorkflowNode } from '../../shared/types';
import { applyInlinePick } from './inlineReferencePick';
import { nodeWithParamValue } from './workflowNodeEdits';
import { WorkflowReferences } from './WorkflowReferences';

// The open inline `{…}` reference menu: which field it hangs under, the
// node/param being edited, and the suggestions to show. Null when idle.
export type InlineReferenceState = {
  anchorEl: HTMLInputElement,
  node: string,
  param: string,
  options: string[],
} | null;

// The inline `{…}` reference autocomplete on a node's parameter field: which
// suggestions to show while typing, what a pick writes into the field, and
// putting the caret back afterwards (the field's value is controlled, so it
// has to be repositioned once React re-renders).
export const useInlineReference = ({
  nodes,
  references,
  onSetNode,
}: {
  nodes: { [name: string]: WorkflowNode },
  references: WorkflowReferences,
  onSetNode: (name: string, node: WorkflowNode) => void,
}): {
  inline: InlineReferenceState,
  closeInline: () => void,
  handleInlineEdit: (nodeName: string, param: string, value: string, caret: number | null, el: HTMLInputElement) => void,
  pickInline: (option: string) => void,
} => {
  const [inline, setInline] = useState<InlineReferenceState>(null);
  const caretRef = useRef<{ el: HTMLInputElement, pos: number } | null>(null);

  const closeInline = () => {
    return setInline(null);
  };

  // Typing / caret moves in a field: refresh the suggestions, or close them when
  // the caret leaves the token.
  const handleInlineEdit = (nodeName: string, param: string, value: string, caret: number | null, el: HTMLInputElement) => {
    const options = references.inlineOptions(nodeName, value, caret);
    if (options === null) {
      setInline(null);
      return;
    }
    setInline({ anchorEl: el, node: nodeName, param, options });
  };

  // Writes a new value into the edited field's parameter and queues the caret to
  // land at `caret` once the controlled input re-renders.
  const commitValue = (nodeName: string, param: string, value: string, caret: number, el: HTMLInputElement) => {
    onSetNode(nodeName, nodeWithParamValue(nodes[nodeName] ?? {}, param, value));
    caretRef.current = { el, pos: caret };
  };

  // Picks the highlighted suggestion: a node leaves the reference open for its
  // output, an output completes it and closes the menu.
  const pickInline = (option: string) => {
    if (inline === null) {
      return;
    }
    const el = inline.anchorEl;
    const picked = applyInlinePick(el.value, el.selectionStart ?? el.value.length, option);
    if (picked === null) {
      setInline(null);
      return;
    }
    commitValue(inline.node, inline.param, picked.value, picked.caret, el);
    if (picked.completed) {
      setInline(null);
    } else {
      setInline({ ...inline, options: references.outputsOf(option) });
    }
  };

  // Restore the caret after a pick rewrote the (controlled) field value.
  useEffect(() => {
    const pending = caretRef.current;
    if (pending !== null) {
      caretRef.current = null;
      pending.el.focus();
      pending.el.setSelectionRange(pending.pos, pending.pos);
    }
  });

  return { inline, closeInline, handleInlineEdit, pickInline };
};
