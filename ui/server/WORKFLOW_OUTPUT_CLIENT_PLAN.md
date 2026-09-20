# Plan (vector 1): client consumes chunks; drop the flat log

Follow-up to `WORKFLOW_OUTPUT_PLAN.md`. That plan kept the client on the old flat
`output` string on purpose. This plan moves the client to the per-task `chunks`
and removes the parent's flat-log rebuild.

## Goal

- Client reads per-task `chunks`, live and when done.
- Server stops building and sending the flat `output` string.
- One shape for output everywhere: an ordered list of `{name, text}` segments.

## Current state

- Server already exposes chunks live: `poll()` returns both `output` and `chunks`
  while running and when done. Chunks arrive live over the result queue.
- Client ignores chunks while running (uses `output`), and only reads chunks when
  the run is done. The parent still rebuilds the flat `output` only for the client.

So this is client-first, then a server cleanup.

## Step 1: client switches to chunks (no server change yet)

- Poller: while running, set the run's chunks from `result.chunks` (not just
  `output`). On done, same. Stop reading `result.output`.
- Store: `WorkflowRun` holds `chunks` as the source of truth. Drop `output` (or
  keep a derived getter that joins chunk text, if a plain-text view is still handy).
  Replace `setRunOutput` with a `setRunChunks` used live.
- Run dialog / button: render from `chunks` in both states (running and done).
  Where a flat view is needed, join the chunk texts.
- Tests: update the vitest specs that assert on `output` to assert on `chunks`.

At this point the client no longer needs `output`, but the server still sends it.

## Step 2: server drops the flat log

- `WorkflowLogBuilder`: remove `_flat` and `output()`; keep the ordered segments and
  `chunks()`.
- `WorkflowRunResult` and `poll()`: drop `output`; return `chunks` (+ `status`,
  `error`).
- `api_models.RunWorkflowResponse`: drop the `output` field.
- The timing line (added as a `WorkflowOutput(BETWEEN, ...)` piece) stays as a final
  chunk, so it still shows.
- Server tests: drop `output` assertions; assert on chunks.

## Interface change

`PollWorkflowResult` loses `output`. This is a breaking change to the poll response,
so ship Step 1 and Step 2 together (or keep `output` one release as deprecated).

## Not in scope

- Live task-status panel UI (#906).
- Any change to how the run is started.

## Checklist

- Follow `ui/client/TEST_UI.md` before reporting the client work done.
