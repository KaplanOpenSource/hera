# Issue #1041 - Decrease view modes (plan)

Source: KaplanOpenSource/hera#1041 (transferred from #1035 task 2). Milestone: Web UI.

**Status: IMPLEMENTED (Option B).** `RawViewToggle.tsx` added, `DocViewSelector.tsx`
deleted, `DetailsViewDocumentContent.tsx` switched to a `rawView` boolean. The plan
below is kept for context.

## The ask (verbatim)

> There really shouldn't be so many view modes for a document. We should have a
> toggle to view raw data that switches between the agent/hermes view and the raw
> document.

> The toggle should be disabled for documents lacking specialized implementations.

So: replace the current multi-button view selector with **one** "Raw" toggle that
flips between a document's specialized/processed view and its raw stored form.

## Current state (what exists today)

- `DocViewSelector.tsx` renders a 4-button `ToggleButtonGroup`: `raw`,
  `formulated`, `agent`, `workflow`. Agent/Workflow buttons are disabled unless the
  doc is of that kind.
- `DetailsViewDocumentContent.tsx` holds the mode in `docView` state
  (`type DocView = 'raw' | 'formulated' | TabKind.Agent | TabKind.Workflow`):
  - `defaultView(doc)`: agent -> `Agent`, workflow -> `Workflow`, else `formulated`.
  - An effect resets `Agent`/`Workflow` back to `formulated` when the doc kind no
    longer matches.
  - Derived flags: `showFormulated = docView !== 'raw'`,
    `showAgentConfig = docView === Agent`, `showWorkflow = docView === Workflow`,
    `showKindEditor = showAgentConfig || showWorkflow`, and `headerHiddenFields`
    (agent hides `resource/type/dataFormat`; workflow hides `type/dataFormat`).
- Separately there is `DetailsVisibilityToggle` (Both / Tree / None). This is a
  DIFFERENT axis (how much of the detail section is shown), NOT a view mode. Out of
  scope - leave it as is.
- Usages of `DocViewSelector` / `DocView` / `docView`: only in
  `DetailsViewDocumentContent.tsx`. No tests depend on `DocViewSelector`
  (`deleteDescField.test.tsx` only mentions "formulated" in a comment).

## Target behavior

One control - a labeled "Raw" toggle (a MUI `Switch`, matching the mock's "Raw View"):

| Doc kind | Toggle | OFF (default)                    | ON (raw)                 |
|----------|--------|----------------------------------|--------------------------|
| Agent    | enabled| agent editor                     | raw document (all fields)|
| Workflow | enabled| workflow/hermes editor           | raw document (all fields)|
| Regular  | see Q1 | formulated field tree            | raw document (all fields)|

The `formulated` mode disappears as a user-facing choice: for agent/workflow docs
the OFF state is their editor; for regular docs the OFF state is the formulated tree.

## Open question (needs a decision before coding)

**Q1. Is the Raw toggle disabled for regular documents? -> DECIDED: Option B.**
The Raw toggle is **always enabled**, for every document kind. A regular doc flips
formulated <-> raw; an agent/workflow doc flips its editor <-> raw. No kind-based
`disabled`.

## Proposed design

Replace the `docView` enum state with a single boolean.

- New state: `const [rawView, setRawView] = useState(false)` (default: specialized).
  - Reset to `false` when the open document changes (replaces the `defaultView`
    effect). No kind-mismatch effect needed - `rawView` is kind-agnostic.
- Derived flags (drop-in replacements for the current ones):
  - `showAgentConfig = !rawView && isAgent`
  - `showWorkflow   = !rawView && isWorkflow`
  - `showFormulated = !rawView`  (raw shows everything; OFF hides boilerplate)
  - `showKindEditor = showAgentConfig || showWorkflow`
  - `headerHiddenFields` - unchanged logic, keyed off the two flags above.
  - This reproduces today's agent/workflow/raw rendering exactly; the only change
    is that `formulated` is no longer separately selectable.
- New component `RawViewToggle.tsx`:
  - A `Switch` + "Raw" label (or `FormControlLabel`), placed where `DocViewSelector`
    was in the control row.
  - `checked={rawView}`, `onChange` toggles.
  - `disabled` per Q1 (Option A: `!(isAgent || isWorkflow)`).
  - Tooltip explaining what it does; when disabled, explain why.

## Steps

1. Decide Q1 (default to Option A).
2. Add `RawViewToggle.tsx`.
3. `DetailsViewDocumentContent.tsx`:
   - Remove `docView` state, `defaultView`, and the kind-mismatch effect.
   - Add `rawView` state + reset-on-doc-change effect.
   - Rewrite the three derived flags as above.
   - Swap `<DocViewSelector .../>` for `<RawViewToggle .../>`.
4. Delete `DocViewSelector.tsx` and the `DocView` type (no other users). Keep
   `TabKind` (used widely elsewhere).
5. Verify `deleteDescField.test.tsx` still passes (it relies on the default,
   non-raw view hiding `HIDE_ON_DESC` fields - still true when `rawView === false`).
6. Consider a small unit test for `RawViewToggle` (enabled/disabled + toggle) and
   maybe assert a regular doc renders formulated by default.
7. Run the full TEST_UI checklist (tsc -> tests -> build).

## Edge cases

- A doc that is both agent-ish and workflow-ish: not possible today
  (`isAgent`/`isWorkflow` are derived from distinct shapes); keep agent priority as
  the current `defaultView` does.
- Switching documents while in raw view: reset to specialized (OFF) on doc change so
  each doc opens in its natural view.
- `DetailsVisibilityToggle` (Both/Tree/None) stays and composes with raw view
  unchanged.

## Out of scope

- `DetailsVisibilityToggle` (separate concern).
- The "Module: <cls>" subtitle / big title icon (that's the leftover of #1035 D7).
