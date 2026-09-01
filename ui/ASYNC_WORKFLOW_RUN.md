# Async Workflow Run - Plan

Two concerns got mixed in the first draft and are now split into two steps:

- **Step 1 - monitor an async process.** Turn the blocking run into an async run
  identified by a token. The client polls for *status* only; the full output comes
  back in one piece when the run finishes. No streaming.
- **Step 2 - stream the output** (later). Add live, incremental output while the
  run is still going (0.5s updates).

This doc plans **step 1** in full and sketches step 2 as a follow-up. Step 1 is the
foundation for issue #906 (per-task progress view).

## Why step 1 first

Doing status-only polling first keeps the design small:

- Output is read once, after the run is done and the reader thread has finished.
  Nothing reads the buffer while it is being written, so **no lock on the output**.
- No byte/line offset, no cursor, no snapshotting. `poll` just reports status and,
  when done, the whole output blob.

Streaming (step 2) is what forces the offset cursor and the concurrent buffer. We
defer that until step 1 works.

## Current state (why this is needed)

- `POST /run_workflow` blocks: the server forks a child, waits on it, and returns
  `{dispatch_id, output}` only after the run finishes.
- Output is captured in memory by `PipeTee` (an OS pipe + reader thread) and
  returned as one string at the end.
- Runs are serialized by a `threading.Lock` in `WorkflowRunner`.
- The Luigi run is a separate subprocess (`python3 -m luigi ... --local-scheduler`).
- On a hard server crash, the forked child (and its Luigi subprocess) are orphaned
  and keep running.

## Step 1 - async run + status polling

### End results

1. `POST /run_workflow` returns immediately with a token; the run happens in the
   background.
2. The client polls status until the run is done, then shows the full output.
3. Subprocesses stop when the server crashes (no orphaned Luigi runs).

### Scope

In scope: async run, token, status polling, output-on-done, busy handling, crash
cleanup. Out of scope: live/streaming output (step 2), concurrent runs, cancellation.

### Design decisions

1. **Mint the token up front.** Today `dispatch_id` is only known after the run
   finishes. We generate our own token at start and return it immediately; the
   workflow's `dispatch_id` is stored into the job later.
2. **Output returned whole, on done.** No incremental delivery. `poll` returns the
   full output only once status is `done`. `PipeTee` stays as it is today
   (`result()` after the reader drains).
3. **Status field.** `poll` reports `running` / `done` / `error`, plus an `error`
   message on failure. Written once by the background thread, read by `poll`; a
   single field assignment, safe under the GIL without a lock.
4. **Busy = reject.** One run at a time. A start request while a run is in progress
   returns busy.
5. **Single job slot.** The next start overwrites the previous finished job. At
   most one stale result ever lingers. No TTL, no kill, no dangling handling.

### Server

1. Job slot in `WorkflowRunner`: `{ token, status, output, error, dispatch_id }`.
   `status` is one of running / done / error.
2. `start(projectName, workflowName) -> token`:
   - Non-blocking try-lock. If held, a run is in progress: return busy.
   - Mint a token, spawn the forked run in a background thread, return the token.
   - The lock is held for the whole run (the background thread), freed when it ends.
3. Crash cleanup:
   - Run the fork in its own process group (`setsid`).
   - In the child, set `PR_SET_PDEATHSIG = SIGKILL` so the kernel kills it when the
     server process dies (hard crashes).
   - Shutdown handler kills tracked process groups on clean exit (SIGTERM / normal).
4. Background thread: run the fork; on completion set `status = done` (store output
   + `dispatch_id`) or `status = error` (store `error` + partial output). The lock
   frees when the run ends.
5. Routes:
   - `POST /run_workflow` -> `{ token }`, or `{ status: "busy" }` if a run is running.
   - `GET /run_workflow/{token}` -> `{ status, output, error }`. Unknown token ->
     `{ status: "not_found" }`.
6. API models: start -> `{ token }` or busy; poll -> `{ status, output, error }`.

### Client

1. `runWorkflow` returns the token instead of the final output.
2. Poll loop (modeled on `ServerReadyGate`'s `setTimeout` + `cancelled` guard): hit
   `GET /run_workflow/{token}` every 500ms, stop when status is done or error.
   Clean up on unmount / dialog close.
3. `done` -> show the full output. `error` -> show the error. `not_found` (server
   restarted or slot overwritten) -> stop polling gracefully. `busy` from start ->
   show a busy message, do not poll.
4. `RunWorkflowButton` `doRun`: start -> get token -> poll status -> show output;
   drive the spinner off poll status instead of a single `await`.

### Tests

Follow `ui/client/TEST_UI.md`. Vitest tests for the poll loop (stop on done/error,
busy, unknown-token stop) and the server job slot (start/poll, busy, not_found,
slot overwrite). Update `ui/client/tests/execPython-coverage.md` if needed.

## Step 2 - stream the output (later)

Add live output while the run is still going, on top of step 1:

- Store output as a list of lines; `poll` gains a `since` line index and returns
  only new lines plus a `next_offset`.
- The reader thread writes lines while `poll` reads them, so the buffer needs a
  small lock (or a snapshot copy). This is the lock step 1 avoids.
- Client appends new lines each poll and advances the offset; on dialog reopen it
  polls from `since=0` to rebuild the full output.
- Progress output that uses carriage returns (tqdm, some Luigi bars) will not
  stream line by line, since a line only lands once its newline arrives. Acceptable.

## Follow-ups (later, not this plan)

- Per-task status + per-task logs via marker injection in the hermes task mixin
  (issue #906).
- Concurrent runs (relax the lock, key all state by token).
- Cancellation (a stop endpoint; process groups already make this easy).
- Busy polling: before a user starts a run, poll to show if the server is busy.
- Hang recovery is out of scope: a hung run holds the slot until server restart.
