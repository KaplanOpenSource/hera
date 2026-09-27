# Plan: one place that defines a reference kind

Adding "reference another node's input" (#1001) costs about 250 lines today. Almost none of it is
the feature. It is the same fact - what a reference looks like - restated in six places. This plan
moves those six behind one class, so the next kind is one small subclass.

Do this refactor first, then add the input kind on top of it.

## The six places

A reference is written `{node.output.key}`. Six bits of code know that shape on their own:

1. `dataflowReference` writes it.
2. The `REFERENCE` regex reads it out of a parameter value.
3. `tokenAtCaret` reads a half-typed one while the user types.
4. `outputHandleId` / `inputHandleId` name the dots a line runs between.
5. `buildDataflowEdges` and `parseDataflowEdgeId` encode it again as an edge id.
6. `WorkflowReferences` and `WorkflowContextMenu` list what a field may point at.

A second kind means editing all six. That is the whole cost.

## The design

Three types, in a new `workflowReferenceKinds.ts`.

### `Reference` - one reference, as a value

Node, kind, key. It replaces the loose `(sourceNode, outputName)` pairs that are passed around now.
It is what flows through the menus and the edit helpers, so no one below the UI handles a section
string by hand.

```ts
export class Reference {
  constructor(readonly node: string, readonly kind: ReferenceKind, readonly key: string) {}

  toString(): string { ... }                                   // {node.section.key}
  edgeIdTo(target: string, param: string): string { ... }
  handleId(): string { ... }
}
```

### `ReferenceKind` - what differs between kinds

An abstract base. Three fields and one method are abstract; everything else is shared, because every
kind formats, anchors and clears the same way.

```ts
export abstract class ReferenceKind {
  abstract readonly section: string;      // the written form
  abstract readonly label: string;        // shown in menus
  abstract readonly handleMark: string;   // the source dot's id fragment

  // The only real per-kind behaviour: what this node offers of this kind.
  abstract namesOf(node: WorkflowNode, catalog: NodeCatalogEntry[]): string[];

  format(node: string, key: string): string { ... }
  handleId(node: string, key: string): string { ... }
  edgeId(node: string, key: string, target: string, param: string): string { ... }
  clearToken(node: string, key: string): RegExp { ... }
  couldBe(typedSection: string): boolean { ... }   // for a half-typed token
  sectionMatch(): string { ... }                   // what the parser accepts
}
```

Two subclasses:

- `OutputReferenceKind` - section `output`. Overrides `namesOf` (the catalog's outputs) and
  `sectionMatch` (it also reads the older `parameters` spelling).
- `InputReferenceKind` - section `Execution.input_parameters`. Overrides `namesOf` only.

`Execution.input_parameters` is the only form Hermes resolves. A node saves its inputs nested under
`Execution`, so `{A.parameters.key}` raises a `KeyError`. No Python change is needed - the dependency
scan already treats `{A.…}` as a requirement on A.

### `ReferenceKindRegistry` - the operations across kinds

Some questions are about all kinds at once, not one kind. They become methods on a class that holds
the list. It builds the combined parse regex once, in its constructor.

```ts
export class ReferenceKindRegistry {
  constructor(private readonly kinds: ReferenceKind[]) { ... }

  all(): ReferenceKind[]
  parseAll(node: string, value: string): Reference[]    // every reference in a parameter value
  ofHandle(handleId: string): Reference | null          // which kind a dot belongs to
  ofEdgeId(edgeId: string): { ref, target, param } | null
  couldBe(typedSection: string): ReferenceKind[]        // kinds a half-typed section still fits
  keysOf(node: string, ...): Reference[]                // everything one node offers
}

export const knownKinds = new ReferenceKindRegistry([new OutputReferenceKind(), new InputReferenceKind()]);
```

Names are kept apart on purpose: `ReferenceKind`, `ReferenceKindRegistry`, `knownKinds`. No type and
value differing only by case or plural.

## What each of the six becomes

| Today | After |
|---|---|
| `dataflowReference` | `ref.toString()` |
| `REFERENCE` regex | `knownKinds.parseAll(...)` |
| `tokenAtCaret` section handling | `knownKinds.couldBe(token.section)` |
| `outputHandleId` / `inputHandleId` | `ref.handleId()`; the target dot stays as it is |
| edge id build + parse | `ref.edgeIdTo(...)` / `knownKinds.ofEdgeId(...)` |
| the two option lists | `knownKinds.keysOf(...)`, one list, labelled by `kind.label` |

## What stays manual

Be honest about the limits.

- The **target** dot (`:in:`) is not per-kind. A reference always lands on a parameter, so only the
  source side varies.
- The **dot itself** is JSX on a row. A kind that needs a new dot in a new place still costs a few
  lines in the row component. The registry cannot render it.

So a third kind is one subclass plus a dot if it needs one. Small, not free.

## Order of work

1. Add `workflowReferenceKinds.ts` with the three types and `OutputReferenceKind` only. Behaviour
   unchanged.
2. Move `workflowDataflow.ts` onto it - writer, parser, handle ids, edge ids. Existing tests should
   pass with only the edge id format changed.
3. Move `WorkflowReferences`, `inlineReferencePick` and `WorkflowContextMenu` onto `Reference`,
   dropping the duplicate option lists.
4. Add `InputReferenceKind`. This should be the subclass and nothing else.
5. Add the right-hand dot on top-level parameter rows in `WorkflowNodeInputs`, and a render slot
   after the value editor in `DetailsViewItem` (threaded through `DetailsViewItemsInObject`).
6. Block a node referencing itself, in `WorkflowCanvasEdits.canConnect` and `buildDataflowEdges`.
7. Tests, then the checklist in `ui/client/TEST_UI.md` in order.

Steps 1-3 are the refactor and show nothing to the user. Steps 4-6 are the feature.

## The trade-off

The refactor is about the same size as doing the feature without it, and it is invisible when it
lands. It pays off on the third kind of reference, not the second. If the input kind is the last one
we will ever add, skip this plan and do the feature directly.

## Open question

An input reference reads a value the other node was *given*, not one it produced, yet Hermes still
makes the reader wait for that node to finish. The link is safe but adds ordering the user may not
expect.
