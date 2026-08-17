# Async Workflow Run - Plan

Goal: change workflow runs from one blocking HTTP request into an async run
identified by a token, with the UI polling for live output. This is the
foundation for issue #906 (per-task progress view), but this plan covers only
the async run itself - not the per-task view and not concurrency.

## End results

1. Subprocesses stop when the server crashes (no orphaned Luigi runs).
2. The output shown in the UI updates every 0.5s while the run is in progress.

## Current state (why this is needed)

- `POST /run_workflow` blocks: server forks a child, waits on it, returns
  `{dispatch_id, output}` only after the run finishes.
- Output is captured in memory by `PipeTee` (an OS pipe + reader thread) and
  returned as one string at the end.
- Runs are serialized by a `threading.Lock` in `WorkflowRunner`.
- The Luigi run is a separate subprocess (`python3 -m luigi ... --local-scheduler`),
  so all output comes back as one shared stdout/stderr stream.
- On a hard server crash, the forked child (and its Luigi subprocess) are
  orphaned and keep running.

## Scope

In scope: async run, token, polling, live output, crash cleanup.
Out of scope: per-task status view (#906), concurrent runs, cancellation.

## Things we also need (not just the two end results)

1. Mint the token up front. Today `dispatch_id` is only known after the run
   finishes. For async we generate our own token at start and return it
   immediately; the workflow's `dispatch_id` is stored into the job later.
2. Incremental output. `PipeTee` currently returns text only at the end (it
   joins the reader thread). For 0.5s updates it needs a snapshot method that
   returns captured output from a byte offset onward, without joining.
3. A done signal. The poll response needs a `status` field
   (running / done / error) so the client knows when to stop polling, plus an
   `error` message on failure.
4. Busy behavior. Decided: reject. There can be only one run at a time; a start
   request while a run is in progress returns busy.
5. Finished-job cleanup. The lock frees when the run process ends, so a new run
   can start right away - no busy-forever, no kill needed. The only leftover is
   one finished job's output + done status the client may not have fetched yet;
   it just sits in the single slot and is overwritten when the next run starts.
   At most one stale output ever lingers. No TTL, no kill, no dangling handling.

### Output delivery: incremental (offset cursor)

Outputs can be long, so we send only the new part each poll instead of the whole
blob:

- Poll: `GET /run_workflow/{token}?since=<offset>`.
- Response: `{ status, chunk, next_offset, error }` where `chunk` is the output
  after `offset`.
- Client appends `chunk` and stores `next_offset` for the next poll.
- On dialog reopen / remount, poll with `since=0` to rebuild the full output.

Single run + single dialog means no ambiguity about the read position, so this
does not add confusion.

## Plan

### Server

1. Job registry in `WorkflowRunner`: `token -> { status, output_so_far, error,
   process, process_group, dispatch_id }`. `status` is one of running / done /
   error.
2. `start(projectName, workflowName) -> token`:
   - Busy is a non-blocking try-lock. If the lock is held, a run is in progress:
     return a busy JSON response. The client shows a "busy" message.
   - Mint a new token.
   - Spawn the forked run in a background thread and return the token
     immediately.
   - The lock is held for the whole run lifetime (the background thread), not
     the request. It frees when the run ends.
3. Crash cleanup:
   - Run the fork in its own process group.
   - In the child, set `PR_SET_PDEATHSIG = SIGKILL` so the kernel kills it when
     the server process dies (covers hard crashes).
   - Add a shutdown handler that kills tracked process groups on clean exit
     (covers SIGTERM / normal shutdown).
4. `PipeTee.read_from(offset)`: return captured output from `offset` onward as
   text (plus the new offset) without joining the reader thread, so output can
   be read while the run is still going.
5. Background thread: run the fork; on completion set `status = done` (store
   `dispatch_id`) or `status = error` (store `error` + partial output). The lock
   frees when the run ends, so the next start can proceed.
6. Routes:
   - `POST /run_workflow` -> `{ token }`, or a busy JSON if a run is in progress.
   - `GET /run_workflow/{token}?since=<offset>` -> `{ status, chunk, next_offset,
     error }`. Unknown token -> not found.
7. Next start overwrites the single slot, discarding any previous finished job.
8. New API models: start -> `{ token }`; poll -> `{ status, chunk, next_offset,
   error }`; plus the busy and not-found shapes.
9. Thread safety: the reader thread writes output while the poll endpoint reads
   it, so guard the buffer/job with a small lock or read a snapshot copy.

### Client

1. `runWorkflow` returns the token instead of the final output.
2. Poll loop (modeled on `ServerReadyGate`'s `setTimeout` + `cancelled` guard):
   hit `GET /run_workflow/{token}?since=<offset>` every 500ms, append `chunk` and
   advance the stored offset each time, stop when `status` is done or error.
   Clean up on unmount / dialog close.
3. Unknown token from poll (server restarted, or slot overwritten): stop polling
   gracefully and show the output collected so far. This also covers a server
   restart mid-run.
4. Busy response from start: show a "busy" message, do not start polling.
5. `RunWorkflowButton` `doRun`: start -> get token -> poll -> show live output;
   drive the spinner off poll status instead of a single `await`.

### Tests

Follow `ui/client/TEST_UI.md`. Add vitest tests for the poll loop (append,
offset advance, stop on done/error, unknown-token stop) and the server job
registry. Update `ui/client/tests/execPython-coverage.md`.

## Follow-ups (later, not this plan)

- Per-task status + per-task logs via marker injection in the hermes task mixin
  (issue #906).
- Concurrent runs (relax the lock, key all state by token).
- Cancellation (a stop endpoint; process groups already make this easy).
- Busy polling: before a user starts a run, poll to show if the server is busy.
- Hang recovery is out of scope: a hung run holds the slot until server restart.
