import { createTheme } from '@mui/material';
import { useMemo } from 'react';
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
