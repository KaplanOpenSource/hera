# Step 2 - Stream workflow output (Option A)

Goal: show the run output while it is still going. Now we only show it at the end.
Builds on Step 1 (start + poll + shared store + poller). See ASYNC_WORKFLOW_RUN.md.

## Approach (Option A: send the whole text each poll)

Each poll returns the full output so far. The client replaces its copy.
No offsets, no line diffing. Fine for normal (small) logs.

Why this fits: `WorkflowLogView` takes one big string and redraws it. The store
keeps the run for the whole app, so closing the dialog loses nothing.

## Server (ui/server)

1. `pipe_tee.py`: add `snapshot() -> str`.
   - Return the text collected so far, without waiting for the run to end.
   - Copy the buffer first: `b"".join(list(self._collected)).decode(errors="replace")`.
   - The `list(...)` copy makes reading safe while the reader thread appends. No lock.
2. `workflow_runner.py`: keep the live buffer on the job so `poll` can read it.
   - `_background` creates the `PipeTee`, stores it on the job (e.g. `job["tee"]`),
     runs the fork, then sets the final output + status.
   - Move the fork/wait code so `_background` holds the tee (refactor `run()` into a
     helper that takes a tee, or inline it). Keep the single-run lock.
   - `poll`: while running, `output = job["tee"].snapshot()`; when done, `output` is
     the final text. Same response shape (`status`, `output`, `error`). `output` is
     just filled during `running` too.
3. Tests: `poll` returns partial output while running; full output when done.

## Client (ui/client)

4. `stores/useWorkflowRunStore.ts`: add an action to set output while still running
   (e.g. `setRunOutput(workflowName, output)`), keeping status = running.
5. `components/workflow/WorkflowRunPoller.tsx`: on every poll (running and done),
   write `output` to the store, not only on done.
6. `components/workflow/log/WorkflowOutputDialog.tsx`: show the growing log while
   running (with a small "running" hint), instead of only a spinner.
7. Tests: poller updates output while running; dialog shows partial output.

## Validate

Follow ui/client/TEST_UI.md: tsc, tests, build. Run server tests in the container.

## Caveat

Progress bars that use `\r` (tqdm, some Luigi bars) show as one long growing line.
Acceptable.

## Later (not now)

Option B (send only the new part with an offset) if logs get large.
