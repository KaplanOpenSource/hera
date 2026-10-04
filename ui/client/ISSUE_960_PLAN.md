# Issue 960 - Tree interactions like a file explorer

Goal: make the Workspace Explorer behave like a file tree.

- Single left click selects a row. It does not open a tab.
- Double left click opens the document tab.
- Ctrl / Shift click keep working (MUI already does this).
- Clicking a branch row selects all documents under it.
- Clicking empty space in the panel clears the selection.
- Right click opens a menu with Delete and Duplicate.
- Double click on a repository asks to pull its contents into the project.

## Where the code is

- `src/components/project/ProjectTreeView.tsx` - owns both trees, selection
  state, URL sync. `handleSelectionChange` currently calls `onSelectItem`, so a
  single click opens a tab. This is the main file to change.
- `src/components/project/DocumentSplitGroup.tsx` - renders branch rows
  (`TreeItem`) and leaf rows (`ProjectDocumentItem`).
- `src/utils/splitTree.ts` - `SplitTree`, `collectBranchKeys`.
- `src/components/project/RepoTreeWhole.tsx` - the repositories tree.
- `src/components/repo/LoadRepositoryButton.tsx` - already has the Python call
  that loads a repository into the project. Reuse it.
- `src/components/project/DeleteSelectedButton.tsx` - already has the bulk
  delete Python call. Reuse the same code shape.

## 1. Split "select" from "open"

In `ProjectTreeView.tsx`:

- `handleSelectionChange` only sets `selectedIds` and clears `repoSelectedIds`.
  Remove the `navigateToItem` and `onSelectItem` calls from it.
- Add `lastClickedRef = useRef<string | undefined>()`. Set it to the newly
  added id on every selection change (both trees).
- Add `openItem(id)` which does `navigateToItem(id)` then `onSelectItem(id)`.
- Put `onDoubleClick` on the documents `SimpleTreeView`. It calls
  `openItem(lastClickedRef.current)`. The browser fires the click first, so the
  selection is already correct and no prop plumbing into the rows is needed.
- Keyboard: Enter also opens. Add `onKeyDown` on the tree, open on `Enter`.

`selectDocument` (used after "Add document") still selects *and* opens. Keep it.

The URL stops following a plain click. That is wanted: the URL names the open
document, not the highlighted one. The two `useEffect`s that sync `selectedIds`
from `docId` stay as they are.

## 2. Branch click selects its documents

In `splitTree.ts` add:

```ts
export const collectSubtreeKeys = (node: SplitTreeNode): string[] => { ... }
```

It returns the branch's own `itemKey`, the keys of nested branches, and
`idDocId(doc.docid)` for every leaf below.

In `ProjectTreeView.tsx`, when the clicked id is a split id
(`isSplitId` from `shared/idDocId`), replace the ids MUI gave us with
`collectSubtreeKeys(thatNode)` merged into the selection. Find the node with a
small `findNode(itemKey)` method on `SplitTree`.

Rule: a plain click on a branch replaces the selection with that subtree.
Ctrl click adds the subtree. Shift click is left to MUI.

With this, the delete button in the project Actions menu works on a whole
branch, which is the main point of the request.

## 3. Background click clears the selection

Wrap the current `<>...</>` in `ProjectTreeView` in a `Box` with
`height: '100%'` and an `onClick`. If
`!(e.target as HTMLElement).closest('.MuiTreeItem-content')` then clear both
`selectedIds` and `repoSelectedIds`. Do not touch the URL.

## 4. Right click menu

New file `src/components/project/TreeContextMenu.tsx`. One component, a MUI
`Menu` anchored at the mouse position (`anchorReference="anchorPosition"`).

- `ProjectTreeView` keeps `contextMenu: { x, y } | null`.
- `onContextMenu` on the documents tree: `preventDefault()`, read the row from
  `closest('.MuiTreeItem-content')`. If the row is not selected, make it the
  only selection first. Then open the menu.
- Items:
  - **Delete** - same confirm + Python as `DeleteSelectedButton`. Pull that
    call out of the button into `deleteDocuments.ts` so both use it.
  - **Duplicate** - new `buildDuplicateDocumentCode.ts`. For each selected
    document, read it and add a copy with `datasourceName` suffixed ` copy`.
    Verify the exact datalayer call before writing it; `All.getDocumentByID`
    plus an `addDocument` on the same collection is the expected shape.
  - Disable both when nothing is selected.

Reload with `fetchProjectDetails(projectName)` after either action, as the
delete button does today.

## 5. Repository double click

In `RepoTreeWhole.tsx`, the row ids are repo ids (`idRepoId`). The parent
already owns repo selection, so:

- Move the repository Python load call from `LoadRepositoryButton.tsx` into
  `loadRepositoryIntoProject.ts` and have the button call it.
- `ProjectTreeView` gets `onDoubleClick` on the repositories tree. Take
  `lastClickedRef.current`, map with `idFromRepoId`. If it is a repo name, ask
  with `useConfirm` ("Load repository X into project?") and on yes call
  `loadRepositoryIntoProject(name)`.
- A double click on a non-repo row there (central folder, a JSON value) does
  nothing.

## 6. Cleanup in the touched file

Remove the three stray `console.log` lines in `ProjectTreeView.tsx`
(`[focus]`, `toolkits`, `project`).

## Tests

New `tests/treeInteractions.test.tsx`:

- single click selects and does not call `onSelectItem`
- double click calls `onSelectItem` once with the clicked document
- branch click selects every document under it
- click on the panel background clears the selection
- right click on an unselected row selects it and opens the menu
- `collectSubtreeKeys` unit tests in `tests/splitTree.test.ts` if it exists,
  else alongside the above

## Order of work

1. Step 6 (logs), step 1 (select vs open) + its tests.
2. Step 2 (branch selection) + tests.
3. Step 3 (background click) + test.
4. Step 4 (context menu, delete then duplicate) + tests.
5. Step 5 (repo double click) + test.

## Open question

Duplicate needs the right datalayer call for "read one document and add a copy
to the same collection". I will check the Python side before step 4 and report
what it is.

## Before reporting done

Run the full `ui/client/TEST_UI.md` checklist in order.
