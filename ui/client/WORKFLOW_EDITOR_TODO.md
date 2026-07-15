# Workflow editor — remaining work

The ReactFlow editor for Hermes workflow JSONs (#898, under #852). This lists
only what's **left**; everything else (node graph, `requires` edges, type
autocomplete from the catalog, name prefill, validation, outputs column, param
source dots, dataflow lines, drawing/removing dataflow references, column
ordering by references) is done.

## Node content from the catalog

- **Default parameter values on prefill.** Picking a type seeds parameter
  **names** with empty values — the catalog (built on `node_lookup.py`) carries
  names only. Read each type's `jsonForm.json` `Execution.input_parameters` to
  seed real defaults.
- **Typed parameter fields.** Replace the generic name/value tree with typed
  inputs (number / bool / path / dropdown) driven by `jsonForm` `GUI.Properties`
  (`App::PropertyPath`, `App::PropertyBool`, …). Required markers already exist.

## Dataflow (input ↔ output references)

Input values like `{C.output.ggg}` / `{C.parameters.ggg}` are parsed and drawn
as lines, referenced nodes order into earlier columns, and lines are editable by
drawing (drag output→input, X/delete to clear).

- **Polish.**
  - Input target handle sits *inside* the node (tree indentation), not on the
    left edge, so the line ends a bit inside the box. (Currently the line is
    lifted above nodes via `zIndex` so it stays visible.)
  - Handle/edge colour is hardcoded `#1976d2` — use the theme palette.
  - Only **top-level string** parameter values are scanned for references
    (nested objects / lists are not).

## Open decisions

- Form fidelity: generic tree (fast) vs schema-driven typed forms.

## Out of scope (separate issues)

- #904 Run workflow from UI · #905 Sync workflow (disk ⇄ DB) · #906 Hermes
  output / logs.
