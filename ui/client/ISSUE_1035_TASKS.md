# Issue #1035 - Document view design feedback

Milestone: Web UI. Six UI refinements to the document view. Source: KaplanOpenSource/hera#1035.

## Tasks

- [ ] **1. Collapsible metadata card**
  Metadata (`version`, `toolkit`, `type`, ...) takes too much space. Put it in a single collapsible card instead of separate view modes.

- [ ] **2. Simplify view modes to one raw-data toggle**
  Replace the multiple view modes with one toggle that switches between the agent/hermes view and the raw document view.

- [x] **3. Make the project tab static**
  Convert the project tab to a fixed element so it can't be accidentally closed or moved.
  Done: `LayoutModel.ts` `TREE_TAB` now has `enableClose/enableDrag/enableRename: false`,
  and its tabset has `enableClose/enableDrop: false`.

- [ ] **4. Fix dark-mode contrast**
  Buttons don't separate clearly in dark mode. Add a dark teal background for contrast.

- [x] **5. Auto Reload toggle**
  Replace the current unclear auto-reload control with an explicit toggle switch labeled "Auto Reload".
  Done: `AutoReloadToggle.tsx` is now a MUI `Switch`; `DashboardHeader.tsx` shows an
  "Auto-reload" label next to it.

- [ ] **6. Move delete button into the document interface**
  Move the document delete button inside the document view, alongside the other controls.

## Color changes (from the mockup vs current code)

The mockup is a custom dark-navy + teal/cyan theme. Our code has no custom palette
(`src/theme.ts` just toggles stock MUI light/dark), so today it renders stock MUI
blue `#1976d2` on flat dark `#121212`. These are the specific color diffs. Most of
them belong under task 4 (dark-mode contrast / dark teal background).

- [x] **C1. App background: flat gray -> dark navy.**
  Done: `theme.ts` dark mode sets `background.default = #0b1220`, `paper = #111a2b`.

- [ ] **C2. Card/panel background + border.**
  Current surfaces use `background.paper` = `#121212` with `boxShadow: 1`
  (e.g. `DetailsViewDocumentContent.tsx:196`). Mockup uses a slightly lighter
  translucent navy card (~`#111827`) with a subtle 1px border instead of a shadow.

- [x] **C3. Primary accent: blue -> teal/cyan.**
  Done: `theme.ts` `primary.main` = `#22d3ee` (dark) / `#0891b2` (light).
  Flows into buttons, graph edges, selection highlights, output dots, etc.

- [x] **C4. Header bar: solid blue -> dark navy.**
  Done: `DashboardHeader.tsx` header theme `background.paper = #0b1220` (was `#1976d2`).

- [x] **C5. Logo color: white -> teal/cyan.**
  Done: `assets/atom.svg` strokes/fills now `#22d3ee`. "Hera UI" text stays white
  (matches mockup).

- [ ] **C6. "ADD EFFECT" button: blue -> cyan.**
  `EffectsListEditor.tsx:83` is `variant="contained"` = blue `#1976d2`. Mockup
  is a bright cyan/teal button with dark text.

- [ ] **C7. Chips get accent colors.**
  "Threshold" / "MaxConcentration" chips are default grey outlined MUI chips
  today (`WorkflowNodeOutputChip.tsx:19`, `ToolkitDetails.tsx:23,28`). Mockup:
  Threshold = blue/cyan tint, MaxConcentration = teal/green tint, on translucent dark.

- [x] **C8. Auto-reload "on" color.**
  Done: the new `Switch` in `AutoReloadToggle.tsx` uses cyan `#22d3ee` when on.

- [x] **C9. Green status text tone.**
  Done: `UserIndicator.tsx` now uses `#4ade80` (and matching `rgba` for the docker suffix).

- [ ] **C10. Sidebar/explorer panel background.**
  The tree panels come from flexlayout-react's stock `dark.css` (`theme.ts:21-34`).
  Mockup panels match the navy theme, so this needs a css override to align.

Note: warning triangles (amber) already match closely (`warning.main #ed6c02` vs
mockup ~`#f59e0b`); CORS red already matches. Left off the list.

## Notes

- Reporter says these are relatively minor and offered to split them into separate issues.
- Figma mockups were attached. The Figma tool also suggested broader redesigns, which are out of scope here.
