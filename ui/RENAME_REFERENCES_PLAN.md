# Rename a node, keep its references (#992)

## The problem

Renaming a node in the workflow canvas changes the node key, the `nodeList`
entry, and any `requires` that name it (`WorkflowEditor.tsx`, `renameNode`).

It does not change the `{oldName.output.key}` and `{oldName.input.key}` tokens
written inside other nodes' input parameters. Those tokens keep pointing at a
node that no longer exists, so every reference to the renamed node breaks
silently.

## What must change

One extra step in `renameNode`: rewrite the node part of every reference token
in the workflow.

Tokens are not only top-level strings. In the repo's own example workflows, a
parameter value is often a dict or a list with the token buried inside. Counted
over all workflow JSON files: 174 tokens at the top level of a value, 246
deeper (up to 5 levels). So the rewrite must recurse.

It is all one JSON tree, so the simplest correct shape is a single recursive
walk over the whole workflow block that rewrites the token in every string it
meets. No need to track which node or which parameter a string came from. The
token shape `{name.section.key}` is specific enough not to hit anything else.

## Plan

1. Add a method to `ReferenceKindRegistry` that takes a string and returns it
   with `oldName` replaced by `newName` in every token of a known kind. The
   registry already owns the `PARSE` regex and `bySection`, so it is the right
   home. Unknown sections stay untouched, same as `parseAll`.

2. Add a small recursive helper (own file, `shared/` or next to the references
   folder) that walks any JSON value - string, list, dict - and returns a new
   value with every string passed through a given rewrite. Pure, no workflow
   knowledge.

3. In `renameNode`, after the existing key / `nodeList` / `requires` work, pass
   the rebuilt `nodes` through the walk with the registry's rewrite. Write the
   result back with `setBlock`.

## Tests

In `ui/client/tests/`:

- rename with an output reference in a plain string value
- rename with an input reference
- token inside a nested dict, and inside a list
- a value with text around the token, e.g. `-{node.output.x}+1e-6`
- a node with a similar name is not touched
- a non-string value (number, bool, null) survives the walk unchanged

## Out of scope

Deleting a node leaves dangling references too. Same bug, separate issue, needs
a product call on what to do with the token. Not part of this change.
