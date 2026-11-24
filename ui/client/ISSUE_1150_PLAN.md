# Issue #1150 - Boolean fields in the Hermes editor

## Problem

- The editor knows only number, string, null, object and array.
- A real `true` / `false` is shown as a string.
- Editing it saves the text `"true"`, not a boolean.

## Plan

1. Add a `boolean` type to `ItemTypesEnum` in `src/components/details/ItemTypeSelector.tsx`.
2. Make `calcItemType` return `boolean` for `true` / `false`.
3. Make `coerceToType` turn a value into a boolean (e.g. `"true"` -> `true`, else `false`).
4. Give the boolean chip its own color.
5. In `src/components/details/DetailsViewItemSingle.tsx`, show a switch with "TRUE" / "FALSE" for booleans, instead of the text box.
6. Add unit tests for the type check, the conversion, and the switch.
7. Run the `TEST_UI.md` checklist.

## Out of scope

- No server or Python change. The type comes from the value itself.
