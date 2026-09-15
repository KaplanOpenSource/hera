# Issue 1068: lists act as dictionaries

## Problem

The details editor has no array type. Arrays are treated as objects.

- `calcItemType` in `ItemTypeSelector.tsx` uses `typeof val === 'object'`, so arrays report as `object`.
- `DetailsViewItem.tsx` sorts child keys as text, so `[10]` comes before `[2]`.
- Every edit spreads the array into a plain object, so `[a, b]` becomes `{"0": a, "1": b}`.
- "Add item" adds a key named `newItem_1` to the array.

The read-only viewer `JsonTreeNode.tsx` is fine. This is editor only.

## Plan

### Step 1: add an array type

Add `array` to `ItemTypesEnum`. Check `Array.isArray` before `typeof`.
Allow switching between object and array. Values keep their order, keys become indices.

### Step 2: render array children correctly

Keep one render path. Only the child list differs.
Objects use `Object.entries`. Arrays use `value.map`, so order is numeric.

### Step 3: edit arrays as arrays

Rebuild an array on every change. Use `map` to replace and `filter` to delete.
Never spread an array into an object.

### Step 4: read-only indices

Pass no `setItemKey` for array children. This hides rename and delete-by-key.

### Step 5: add and remove items

"Add item" appends to the end of the array.
Delete removes the element and shifts the rest.

## Follow ups (not part of this plan)

- Drag and drop to reorder list elements.
- A button to convert `{"0": a, "1": b}` back into a list. Manual only, since some dicts use numbers as real keys.
