# Issue 1121 - Minor visual changes with the UI

High level plan. Five small changes, all in the header and the dock tab strip.

## 1. Swap auto-reload and dark mode switches

Dark mode is played with often, auto-reload almost never.

- Move the "Dark mode" switch out of the settings dialog into the header,
  where the auto-reload toggle sits now.
- Move the auto-reload toggle into the settings dialog, next to the reload
  interval slider it controls.
- Files: `src/components/header/DashboardHeader.tsx`,
  `src/components/project/ProjectViewSettingsButton.tsx`,
  `src/components/header/AutoReloadToggle.tsx`,
  `src/elements/ThemeModeSwitch.tsx`.

## 2. Make the tab island fullscreen button clearer

The button is flexlayout's built-in maximize button, so it is styled by CSS,
not by our components.

- Add rules for `.flexlayout__tab_toolbar_button-max` /
  `-min` to the dark mode override block in `src/theme.ts`: stronger icon
  color, visible on hover, and always shown instead of only on hover.
- Add the same rules for light mode (the override block is empty there today).
- Keep the position as it is.
- File: `src/theme.ts`.

## 3. Add-project button next to the project droplist

- Move `<AddProjectButton />` from the right hand stack to the left hand
  stack, right after `<ProjectChooser />`.
- File: `src/components/header/DashboardHeader.tsx`.

## 4. Build info and CORS tag to the bottom of the settings menu

- Remove `<VersionShower />` and `<CorsIndicator />` from the header.
- Render them at the end of the settings dialog content, under a divider.
- Files: `src/components/header/DashboardHeader.tsx`,
  `src/components/project/ProjectViewSettingsButton.tsx`,
  `src/components/header/VersionShower.tsx`,
  `src/components/header/CorsIndicator.tsx`.

## 5. Username in the projects droplist header

- The username is fetched inside `UserIndicator`. Pull that fetch into a small
  shared hook or store so the chooser can read it too.
- Show `"<username>'s projects"` as the droplist label, and drop the separate
  username text from the header.
- Keep the "(docker)" mark somewhere, likely in the droplist tooltip.
- Files: `src/components/header/UserIndicator.tsx`,
  `src/components/header/ProjectChooser.tsx`,
  `src/components/header/DashboardHeader.tsx`, plus a new
  `src/stores/useServerUserStore.ts` (or hook) for the fetch.

## Order of work

1. Header moves (items 3, 1, 4) - all touch `DashboardHeader.tsx`, do together.
2. Settings dialog additions (items 1, 4).
3. Username store and droplist label (item 5).
4. Fullscreen button CSS (item 2).

## Checks

Tests that assert header contents may break (version text, username, add
button). Run the full `TEST_UI.md` checklist at the end.
