import { CssBaseline, ThemeProvider, createTheme } from '@mui/material';
import { useMemo } from 'react';
import { BrowserRouter, Routes, Route } from 'react-router-dom';
import { SnackbarProvider } from 'notistack';
import { Dashboard } from './Dashboard';
import { useViewSettingsStore } from './stores/useViewSettingsStore';

export default function App() {
  const themeMode = useViewSettingsStore((s) => s.viewSettings.themeMode);
  const theme = useMemo(() => createTheme({ palette: { mode: themeMode } }), [themeMode]);
  return (
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <SnackbarProvider maxSnack={6} anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}>
        <BrowserRouter>
          <Routes>
            <Route path="/:projectName/:docId" element={<Dashboard />} />
            <Route path="/:projectName" element={<Dashboard />} />
            <Route path="/" element={<Dashboard />} />
          </Routes>
        </BrowserRouter>
      </SnackbarProvider>
    </ThemeProvider>
  );
}
