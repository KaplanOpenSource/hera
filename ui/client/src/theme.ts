import { createTheme } from '@mui/material';
import { useEffect, useMemo } from 'react';
import lightFlexlayoutUrl from 'flexlayout-react/style/light.css?url';
import darkFlexlayoutUrl from 'flexlayout-react/style/dark.css?url';
import { ThemeMode, useViewSettingsStore } from './stores/useViewSettingsStore';

export const buildAppTheme = (mode: ThemeMode) => {
  const isDark = mode === ThemeMode.Dark;
  return createTheme({
    palette: {
      mode,
      // Cyan accent: bright in dark mode, a deeper teal in light mode for contrast.
      primary: { main: isDark ? '#22d3ee' : '#0891b2' },
      // Dark mode gets a navy background instead of MUI's flat #121212.
      ...(isDark ? { background: { default: '#0b1220', paper: '#111a2b' } } : {}),
    },
  });
};

// The global app theme, following the mode chosen in Settings. Dialogs that open from
// the dark header pass this so they follow the app mode instead of the header's theme.
export const useAppTheme = () => {
  const themeMode = useViewSettingsStore((s) => s.viewSettings.themeMode);
  return useMemo(() => buildAppTheme(themeMode), [themeMode]);
};

// flexlayout ships full light/dark stylesheets that both target the same global
// selectors, so only one can be active at a time. Keep a single <link> in the page
// head and point it at the stylesheet matching the app mode.
export const useFlexlayoutTheme = () => {
  const themeMode = useViewSettingsStore((s) => s.viewSettings.themeMode);
  useEffect(() => {
    const id = 'flexlayout-theme';
    let link = document.getElementById(id) as HTMLLinkElement | null;
    if (!link) {
      link = document.createElement('link');
      link.id = id;
      link.rel = 'stylesheet';
      document.head.appendChild(link);
    }
    link.href = themeMode === ThemeMode.Dark ? darkFlexlayoutUrl : lightFlexlayoutUrl;

    // In dark mode, retint flexlayout's panels (which default to black) to the
    // app's navy so the tree/explorer sidebar matches the content area.
    const overrideId = 'flexlayout-theme-override';
    let override = document.getElementById(overrideId) as HTMLStyleElement | null;
    if (!override) {
      override = document.createElement('style');
      override.id = overrideId;
      document.head.appendChild(override);
    }
    override.textContent = themeMode === ThemeMode.Dark
      ? `.flexlayout__layout {
          --color-background: #0b1220;
          --color-base: #0b1220;
          --color-1: #111a2b;
          --color-2: #16233b;
          --color-tab-selected-background: #1f2d47;
        }`
      : '';
  }, [themeMode]);
};
