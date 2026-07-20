import { createTheme } from '@mui/material';
import { useEffect, useMemo } from 'react';
import lightFlexlayoutUrl from 'flexlayout-react/style/light.css?url';
import darkFlexlayoutUrl from 'flexlayout-react/style/dark.css?url';
import { ThemeMode, useViewSettingsStore } from './stores/useViewSettingsStore';

export const buildAppTheme = (mode: ThemeMode) => {
  return createTheme({ palette: { mode } });
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
  }, [themeMode]);
};
