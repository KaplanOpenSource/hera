# Run workflow button

Goal: let the play button optionally save the workflow before running, without
adding more buttons. Use the button's tooltip + a right-click menu.

## Behavior

- Left click: run the workflow (current behavior).
- Tooltip: "Run workflow (right click for more options)".
- Right click: open a small context menu with these items:
  - Run
  - Run with save (save this doc, then run - one time)
  - Always save before run (sets the flag)
  - Stop saving before run (clears the flag; only shown when flag is on)

## The flag

- Stored in local storage, e.g. key `workflow.alwaysSaveBeforeRun` = `"true"`.
- When ON:
  - Left click saves the doc first, then runs.
  - The icon changes from plain play to play-with-save (small save badge/mini icon).
- When OFF:
  - Left click just runs.
  - Plain play icon.

## Save integration

- "Save" = the same action as the header's Done button:
  `setDoc(new DocumentObj(shownDoc, doc.project))`.
- The button needs access to that save action, so pass a `save` callback (and
  an `isChanged` flag) into `RunWorkflowButton`.
- When saving is ON, unsaved changes are fine: the click saves first, so the
  button stays enabled.
- When saving is OFF and there are unsaved changes, a plain run would execute the
  stale saved version, so the button is disabled (reason "Save changes before
  running"). The right-click menu still opens (it lives on a wrapper Box, not the
  disabled button) so the user can pick "Run with save" or turn saving back on.
- The plain "Run" menu item is disabled while there are unsaved changes;
  "Run with save" is disabled when there is nothing to save.

## Files

- `RunWorkflowButton.tsx` - add right-click menu, flag read/write, icon swap,
  save-then-run flow.
- `DetailsViewDocumentContent.tsx` - pass `save` callback + `isChanged` to the
  button (two call sites).

## Notes

- One button only. No extra buttons anywhere.
- Icon swap is the only visual signal that the flag is on.
