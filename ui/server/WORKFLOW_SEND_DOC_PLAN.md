# Plan (vector 2): send the workflow doc from the client

Today the client sends only `projectName` + `workflowName`. The server looks the
document up in the DB (`getWorkflowListDocumentFromDB`). This plan sends the workflow
content from the client so the run skips that DB lookup.

## Goal

- Client posts the workflow JSON (plus its name) when starting a run.
- Server builds the run straight from that JSON, no DB fetch.

## Why it is possible

The build only needs three things from a doc: the workflow JSON
(`desc['workflow']`), the name (`desc['workflowName']`), and a resource path
(`doc['resource']`). None require a live DB document. `prepareWorkflowRunFromDoc`
uses only those. So a synthesized doc works.

## Server changes

- `RunWorkflowPayload`: replace `workflowName` with a `doc` field (the whole
  workflow document). Keep `projectName`.
- The child: instead of `getWorkflowListDocumentFromDB`, wrap the sent doc in
  `SentWorkflowDoc` and feed it to `prepareWorkflowRunFromDoc`. The wrapper exposes
  the dict the way a saved doc is read: `.desc` (attribute) and `['resource']` (item).
- The doc carries its own `resource`, so the run writes the JSON to the same path a
  saved doc would, custom export paths included.

## Client changes

- The editor already holds the whole document. Send it as `doc`, with
  `desc.workflowName` filled from the resolved name first (the doc's desc may leave
  it unset).

## Decisions

- Send the whole doc, not hand-picked fields. Simpler, and it mirrors a real saved
  document. No DB-lookup fallback, so drop `getWorkflowListDocumentFromDB` from this
  path.
- No save needed. Saving to the DB is handled on the UI, separately. The run works
  straight from the sent doc.

## Open questions

- Resolved: the run uses the doc's own `resource` field (what the DB doc carries),
  so a sent-doc run and a saved-doc run of the same workflow write to the same file.

## Not in scope

- Output transport (that is vector 1).

## Status

Done. The client sends the whole `doc` in the start payload; the server wraps it in
`SentWorkflowDoc` and feeds it to `prepareWorkflowRunFromDoc`, with no DB lookup.
