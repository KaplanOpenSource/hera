# Plan: move workflow output onto the result queue

## Goal

Drop `PipeTee`. Send workflow output to the parent as per-task messages on the
result queue, live. Groups come through directly, no boundary markers.

## Why

- Only the in-process path is left. The router always runs.
- One channel for everything instead of a pipe plus a queue.
- Per-task chunks arrive live, not only at the end.

## Current flow

```
child: fd 1/2 -> router capture pipe -> router thread
                    -> per-task buckets (chunks)
                    -> forward bytes -> tee pipe -> parent tee thread -> live text
child: end -> result queue -> parent (one final message)
```

- `PipeTee`: parent-side pipe reader. Live text for `poll()`.
- Result queue: one message at the end (success or error).

## Proposed flow

```
child: fd 1/2 -> router capture pipe -> router thread
                    -> result queue: {task, text} messages, live
child: end -> result queue: done message
parent: queue reader thread -> builds per-task chunks + flat log for poll()
```

- Keep the router capture pipe. It must exist to catch fd 1/2 (incl. outside programs).
- Router tags each chunk with the current task and puts a `{task, text}` message on
  the queue. One channel, so order and grouping stay exact.
- Parent gets a queue reader thread (replaces the tee thread). From the messages it
  builds the per-task chunks AND rebuilds the flat log string in the shape the client
  parses today, so the client does not change now.
- `poll()` returns both, same as today.
- Remove `PipeTee`.

## Client: leave for now

- Keep the current client parsing untouched. The parent hands it the same flat log.
- Later: change the client to consume the per-task chunks live and drop the flat log
  rebuild on the parent. The chunks already arrive live, so this is a client-only
  follow-up. Tracked here, not done now.

## Messages

NamedTuples, a `Union` narrowed by `isinstance`:

```
WorkflowOutput(task, text)                          # live, many
WorkflowDone(dispatch_id, exec_seconds)             # once
WorkflowError(error)                                # on failure
```

## Not handled on purpose

- Hard child crash with no message sent. We accept it. The child's uncaught output
  lands on the runner's stderr. No crash-recovery logic.

## Not in scope

- Live task-status panel UI (#906). This is transport only.
