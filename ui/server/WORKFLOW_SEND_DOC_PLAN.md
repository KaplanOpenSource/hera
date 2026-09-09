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

- `RunWorkflowPayload`: add a `workflow` field (the workflow JSON). Keep
  `workflowName`; keep `projectName`.
- The child: instead of `getWorkflowListDocumentFromDB`, build a doc-like object from
  the payload: `desc = {'workflow': <json>, 'workflowName': <name>}` and a `resource`
  path under the toolkit `FilesDirectory`. Feed that to `prepareWorkflowRunFromDoc`.
- Pick the `resource` path deliberately (where the workflow JSON is written). Match
  what the DB path used so files land in the same place.

## Client changes

- The editor already holds the workflow. Serialize it to the hermes workflow JSON
  shape and send it in the start payload.

## Decisions

- Require the workflow JSON. The client always sends it. No DB-lookup fallback, so
  drop `getWorkflowListDocumentFromDB` from this path.
- No save needed. Saving to the DB is handled on the UI, separately. The run works
  straight from the sent JSON.

## Open questions

- Exact resource path + name rules, so a sent-doc run and a saved-doc run do not
  collide on disk.

## Not in scope

- Output transport (that is vector 1).

## Status

Parked. Revisit after vector 1.
