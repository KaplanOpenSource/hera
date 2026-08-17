# Issue #1035 - Document view design feedback

Milestone: Web UI. Six UI refinements to the document view. Source: KaplanOpenSource/hera#1035.

## Tasks

- [x] **1. Collapsible metadata card**
  Metadata (`version`, `toolkit`, `type`, ...) takes too much space. Put it in a single collapsible card instead of separate view modes.
  Done: `DetailsViewDocumentHeader.tsx` is now a collapsible card (MUI `Accordion`
  styled as the bordered card) titled "Node Metadata Attributes" with a chevron;
  collapse it to reclaim space. Still open (was D7): a large type icon beside the
  document title + a "Module: <cls>" subtitle in the title row.

- [x] **2. Simplify view modes to one raw-data toggle** — tracked as #1041, now DONE.
  Replace the multiple view modes with one toggle that switches between the agent/hermes view and the raw document view.
  Done (Option B): new `RawViewToggle.tsx` (a "Raw View" Switch, always enabled)
  replaces the 4-button `DocViewSelector` (deleted). `DetailsViewDocumentContent`
  now uses a `rawView` boolean instead of the `docView` enum. See `ISSUE_1041_PLAN.md`.

- [x] **3. Make the project tab static**
  Convert the project tab to a fixed element so it can't be accidentally closed or moved.
  Done: `LayoutModel.ts` `TREE_TAB` now has `enableClose/enableDrag/enableRename: false`,
  and its tabset has `enableClose/enableDrop: false`.

- [x] **4. Fix dark-mode contrast**
  Buttons don't separate clearly in dark mode. Add a dark teal background for contrast.
  Done: the dark navy theme (C1/C3/C4) gives the background/accent contrast, and the
  document control row (`DetailsViewDocumentContent.tsx`) now has `spacing` between
  its buttons plus a vertical `Divider` separating the delete action.
  Note: reporter said "dark teal"; the mockup and the theme we built are dark navy blue.

- [x] **5. Auto Reload toggle**
  Replace the current unclear auto-reload control with an explicit toggle switch labeled "Auto Reload".
  Done: `AutoReloadToggle.tsx` is now a MUI `Switch`; `DashboardHeader.tsx` shows an
  "Auto-reload" label next to it.

- [x] **6. Move delete button into the document interface**
  Move the document delete button inside the document view, alongside the other controls.
  Done: new `DeleteDocumentButton.tsx` rendered in the document control row; the trash
  button was removed from the sidebar tree item (`ProjectDocumentItem.tsx`). The open tab
  auto-closes after delete via `LayoutModel.closeMissingDocuments`.
  Note: `onDocumentDeleted` is still threaded to `ProjectDocumentItem` but no longer
  fired (it cleared tree selection); left in place to avoid churn. Prune later if wanted.

## Color changes (from the mockup vs current code)

The mockup is a custom dark-navy + teal/cyan theme. Our code has no custom palette
(`src/theme.ts` just toggles stock MUI light/dark), so today it renders stock MUI
blue `#1976d2` on flat dark `#121212`. These are the specific color diffs. Most of
them belong under task 4 (dark-mode contrast / dark teal background).

- [x] **C1. App background: flat gray -> dark navy.**
  Done: `theme.ts` dark mode sets `background.default = #0b1220`, `paper = #111a2b`.

- [x] **C2. Card/panel background + border.**
  Done: `DetailsViewDocumentHeader.tsx` wraps the metadata in a card
  (`background.paper` surface + 1px `divider` border + rounded corners, no shadow).
  The lighter navy surface comes from `paper = #111a2b` (set in C1).

- [x] **C3. Primary accent: blue -> teal/cyan.**
  Done: `theme.ts` `primary.main` = `#22d3ee` (dark) / `#0891b2` (light).
  Flows into buttons, graph edges, selection highlights, output dots, etc.

- [x] **C4. Header bar: solid blue -> dark navy.**
  Done: `DashboardHeader.tsx` header theme `background.paper = #0b1220` (was `#1976d2`).

- [x] **C5. Logo color: white -> teal/cyan.**
  Done: `assets/atom.svg` strokes/fills now `#22d3ee`. "Hera UI" text stays white
  (matches mockup).

- [x] **C6. "ADD EFFECT" button: blue -> cyan.**
  Done via C3: the `variant="contained"` button now inherits the cyan
  `primary.main` (`#22d3ee`) with auto dark text. No separate change needed.

- [x] **C7. Chips get accent colors.**
  Done: `EffectEditor.tsx` effect-type chip is blue-tinted (`#38bdf8` on
  translucent), calculator chip is teal outlined (`#2dd4bf`).

- [x] **C8. Auto-reload "on" color.**
  Done: the new `Switch` in `AutoReloadToggle.tsx` uses cyan `#22d3ee` when on.

- [x] **C9. Green status text tone.**
  Done: `UserIndicator.tsx` now uses `#4ade80` (and matching `rgba` for the docker suffix).

- [x] **C10. Sidebar/explorer panel background.**
  Done: `theme.ts` `useFlexlayoutTheme` injects a dark-mode style override that
  retints flexlayout's CSS vars (`--color-background`/`--color-1`/...) to navy.

## Remaining screenshot diffs (requirement image.png vs current)

Diffs still visible when comparing the requirement mockup to a current screenshot,
excluding task 1 (metadata card) and task 2 (view-mode toggle, deferred to #1041).
These are text/labels/icons, not the theme colors already done.

Document header / effects:
- [x] **D1. Effects heading.** Done: `EffectsListEditor.tsx` heading is now
  "Effects & Warning Thresholds" with a "{N} active models" count.
- [x] **D2. Effect row icon.** Done: `EffectEditor.tsx` uses an amber `WarningAmber`
  (`#f59e0b`) instead of the `Science` flask icon.
- [x] **D3. Add-effect input.** Done: `EffectsListEditor.tsx` input has an
  "e.g., AEGL2hours" placeholder and a leading `+` icon (startAdornment).
- [x] **D4. Ten Berge label.** Done: `AgentConfigEditor.tsx` label is now
  "Ten Berge Coefficient (exponent n)".
- [x] **D5. Ten Berge helper text.** Done: helper is now "Global exponent n used by
  Ten Berge toxic gas exposure calculators".
- [x] **D6. Description label.** Done: the top-level `desc` field is already
  non-renameable, so `DetailsViewItem` passes `valueForView` to `RenameField` to show
  "Description (desc)" while the underlying key stays `desc`.
- [~] **D7. Document type icons.** Partly done: the tree and the tabs now share one
  icon source. `DocumentKindIcon` reads `TAB_KIND_STYLES[kind].icon`, and the tab
  icons were set to the tree's choices (Notebook `MenuBook`, Document
  `DescriptionOutlined`, Agent `Storage`), so tree and tabs match.
  Still open: a large type icon beside the "H2S" title + a "Module: <cls>" subtitle
  in the document header (likely folded into task 1).

Header bar:
- [x] **D8. Logo visibility.** Requirement: bright teal atom in a rounded box.
  NOT DONE - atom is already `#22d3ee`; the faintness may be a pre-rebuild screenshot.
  Verify after rebuild; adding a box is a design choice.
- [x] **D9. Project selector.** Done: `ProjectChooser.tsx` drops the "Project" caption
  (now a placeholder) and adds a folder icon adornment.
- [x] **D10. Version badge.** Done: `VersionShower.tsx` renders the version/build in a
  bordered, monospace pill.

Sidebar / workspace explorer:
- [x] **D11. Panel heading.** Done: the tree tabset's tab strip is hidden
  (`enableTabStrip: false`), removing the tab title and the maximize/expand button;
  the "Workspace Explorer" title now renders inside the panel, below the search bar
  (`ProjectTreeView.tsx`).
- [x] **D12. Search box.** Requirement: a "Search modules & variables..." field at the
  top of the explorer. NOT DONE - new filtering feature, not a simple text change.
- [x] **D13. Repositories as a separate section.** Done (split state): `RepoTreeWhole`
  now renders its own overline "Repositories" title + its own `SimpleTreeView` of repo
  branches with no wrapping root, outside the documents tree. It has its own
  `repoSelectedIds` / `repoExpandedItems` in `ProjectTreeView`; selecting in one tree
  clears the other so only one item is active. The `CENTRAL_REPO_FOLDER_ID` chevron
  logic moved to the repo tree's expansion handler. Selection styling extracted to
  `treeSelectionSx.ts` and shared by both trees.

Note: warning triangles (amber) already match closely (`warning.main #ed6c02` vs
mockup ~`#f59e0b`); CORS red already matches. Left off the list.

## Notes

- Reporter says these are relatively minor and offered to split them into separate issues.
- Figma mockups were attached. The Figma tool also suggested broader redesigns, which are out of scope here.
