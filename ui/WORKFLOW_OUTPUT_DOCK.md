# Workflow output as a dock panel (design note)

Idea: replace the modal output dialog (`WorkflowOutputDialog`) with a
flexlayout panel so the run output can float, be dismissed with Esc, or be
docked next to the graph for later work.

Status: **not built** — this note is for deciding whether to do it now or later.

## Why

- Modal blocks the canvas; a dock/float panel lets you watch the log while
  editing or running other things.
- flexlayout-react is already the app's layout engine, so this fits the grain.

## What it touches

1. **`LayoutPanel.tsx`** — add `LayoutComponent.WorkflowOutput` and a case that
   renders the log panel from store state.
2. **`LayoutModel.ts`** — add `openWorkflowOutput(...)` that adds or focuses the
   tab (as a floating tab; see below).
3. **New store** (e.g. `useWorkflowRunStore`) holding run state
   `{ running, output, error }`. A dock tab's content is static config, so the
   live run state can't live in the tab — it must live in a store the panel
   reads and `RunWorkflowButton` writes.
4. **`RunWorkflowButton.tsx`** — instead of local dialog state, write run state
   to the store and call the layout action to open/focus the output tab.
5. **Wiring** — the button needs the layout action. The model lives in
   `ProjectLayout`; expose the "open output" action via the store or the context
   that already holds the model.

The existing `WorkflowLogView` (classified log + filters + copy) is reused
unchanged as the panel body. `WorkflowOutputDialog` would be removed.

## Floating + Esc + docking

- flexlayout supports floating tabs. Set `tabEnableFloat: true` (global or per
  tab) so the tab can pop out to a floating window; the user can drag it back to
  dock it. This gives "float now, dock for later work" for free.
- **Esc to close**: not built in. Add a `keydown` listener while the output tab
  exists that calls `Actions.deleteTab(outputTabId)` on Escape. Small, but it is
  custom code.

## Open decision: one tab or one per workflow

- **One shared tab** — a single "Workflow output" tab, refocused each run.
  Simpler, less clutter. Store keyed by nothing (just the latest run).
- **One per workflow** — tab id per workflow name, so runs can sit side by side.
  Store keyed by workflow name. More code, more UI.

Recommendation: start with **one shared tab** (floating by default); add
per-workflow later only if needed.

## Rough effort

Medium. ~1 new store, ~1 new panel case, ~2 edits (LayoutModel,
RunWorkflowButton), plus the Esc handler and float config. No server changes.
The current dialog stays working until the swap is done.

## Fallback

If dock/float turns out fiddly, a lighter middle ground is a non-modal MUI
`Drawer` (bottom) instead of a `Dialog` — keeps the canvas usable without
touching the layout engine. Not as flexible as a real dock tab.
