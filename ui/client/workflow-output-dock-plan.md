# Workflow run output as a dock tab

## Why

A workflow run's output shows in a modal dialog. The modal blocks the rest of the UI.
Its open state lives inside the run button, and there are two run buttons for the same
workflow, so only the clicked one opens.

Move the output into a normal dock tab, like the workflow canvas already is. Then it can
sit beside the canvas, be resized and moved, and stay open while you work.

## Decisions

- The dialog is replaced, not kept alongside.
- The output tab docks at the bottom of the details panel, in its own tabset.
- The tab opens on its own when a run starts.

## Approach

The run store already holds every run, keyed by workflow name, and the poller already
fills it. The layout watches that store and opens the tab. The run button then owns no
"is the output showing" state at all.

## Steps

1. **Split the dialog.** Pull its body into a plain view component. The modal shell goes
   away; the running hint, error text, log and filter toolbar stay as they are.
2. **Add an output panel.** It takes a workflow name, reads that run from the store, and
   renders the view. No run yet means a short "No run yet." line.
3. **Add the tab to the layout model.** A new tab kind with its own id prefix, following
   the canvas tab: first one docks under the details panel, later ones join that tabset.
4. **Open it from the layout.** Watch the run store; when a workflow starts running, open
   or focus its output tab.
5. **Strip the dialog from the run button.** Start failures keep using the snackbar,
   since a failed start makes no run and so has no tab.
6. **Tests.** One test file for the tab behaviour, mirroring the canvas tab tests. One for
   the panel: live output, running hint, error state. Update any test that asserts on the
   dialog.

## Areas touched

- Workflow log components (`src/components/workflow/log/`)
- Dock layout (`src/components/layout/`)
- The run button (`src/components/workflow/RunWorkflowButton.tsx`)
- Tests (`tests/`)

## Verification

Run the full `ui/client/TEST_UI.md` checklist in order: type-check, unit tests, build.

Then check by hand:

- Run a workflow. An output tab appears under the details panel and fills live.
- Run a second one. It joins the same tabset; the details panel keeps its size.
- Close the tab mid-run, run again. The run is unaffected and the tab comes back.
- Run from both the canvas and the details toolbar. Both reach the same tab.
- Make a run fail. The error shows above the log, and the log stays.
- Drag the tab elsewhere. It still updates.
