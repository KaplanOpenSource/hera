# Low-level steps: rename a node, keep its references (#992)

Scope: the plan in `RENAME_REFERENCES_PLAN.md`. Three steps. Each one builds and
passes tests on its own.

## Step 1 - `utils/mapStrings.ts`

New file, one exported function, nothing imports it yet.

```ts
export const mapStrings = (value: any, rewrite: (text: string) => string): any => {
```

- A string returns `rewrite(value)`.
- An array returns `value.map(item => mapStrings(item, rewrite))`.
- A plain object returns a new object with the same keys, each value mapped.
  Keys are not rewritten.
- Anything else (number, boolean, null, undefined) returns `value` as is.
- Use `Array.isArray` first, then `typeof value === 'object' && value !== null`.
  No `reduce` - a `for...of` over `Object.keys` with an accumulator object.

New test file `tests/mapStrings.test.ts`:

- a bare string is rewritten
- a nested dict and a nested list are rewritten at depth
- keys are untouched even when a key matches
- numbers, booleans, `null` survive unchanged
- the input object is not mutated (compare against a frozen copy)

## Step 2 - `ReferenceKindRegistry.renamedNode`

The registry already owns `PARSE` and `bySection`, so the token rewrite belongs
there.

```ts
// One value's text with every known reference to `oldName` pointing at `newName`.
renamedNode(value: string, oldName: string, newName: string): string {
```

- `value.replace(PARSE, (match, node, section, key) => { ... })`.
- Return `match` unchanged when `node !== oldName`, or when
  `this.bySection(section)` is null. An unknown section stays as written, same
  as `parseAll` skips it.
- Otherwise return `` `{${newName}.${section}.${key}}` ``.
- Keep the section text exactly as it was matched. That way the older
  `{a.parameters.x}` spelling stays that spelling and only the name moves.
- Whitespace inside the braces is dropped (`{ a.output.x }` becomes
  `{x.output.x}`). Acceptable, and worth a test so it is a decision, not a
  surprise.
- `String.replace` with a `/g` regex resets `lastIndex` itself, so sharing
  `PARSE` with `parseAll` is safe.

Extend `tests/ReferenceKindRegistry.test.ts`:

- an output reference is renamed
- an input reference (`{a.Execution.input_parameters.k}`) is renamed
- the older `parameters` spelling keeps its spelling
- a different node name in the same string is left alone
- a name that only shares a prefix (`ab` when renaming `a`) is left alone
- two tokens in one string are both renamed
- an unknown section is left alone
- text around the token survives: `-{a.output.x}+1e-6`

## Step 3 - wire it into `renameNode`

`WorkflowEditor.tsx`, in `renameNode`, after the existing loop that rebuilds
`nodes` with the new key and the renamed `requires`:

- `const renamed = mapStrings(nodes, text => knownKinds.renamedNode(text, oldName, name));`
- pass `renamed` to `setBlock` instead of `nodes`.

Put the walk after the key loop, not inside it, so it sees the whole node map in
one pass.

`requires` is a plain string, not a token, so the walk cannot touch it - the
loop has already fixed it and `renamedNode` leaves bare names alone. Worth one
test that pins this.

Extend `describe('WorkflowEditor renameNode')` in `tests/WorkflowEditor.test.tsx`:

- a reference in a top-level string parameter is renamed
- a reference nested in a dict parameter is renamed
- a reference inside a list parameter is renamed
- a `requires` naming a third node is still untouched
- the renamed node's own parameters are rewritten too (a self reference, and a
  reference it holds to another node, which must not change)

## Then

Run the full `ui/client/TEST_UI.md` checklist before calling it done.

## Out of scope

Deleting a node leaves dangling references. Separate issue, needs a product call
on what to do with the leftover token.
