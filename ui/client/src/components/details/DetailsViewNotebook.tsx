import { Box, CircularProgress, Typography } from '@mui/material';
import { useEffect, useState } from 'react';
import { serverApi } from '../../api/api';
import { ThemeMode, useViewSettingsStore } from '../../stores/useViewSettingsStore';

export const DetailsViewNotebook = ({
  rootDir,
  resource,
}: {
  rootDir: string,
  resource: string,
}) => {
  const [notebookUrl, setNotebookUrl] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const themeMode = useViewSettingsStore((s) => s.viewSettings.themeMode);

  useEffect(() => {
    let cancelled = false;
    setNotebookUrl(null);
    setError(null);

    (async () => {
      const data = await serverApi.ensureJupyter(rootDir, themeMode === ThemeMode.Dark);
      if (cancelled) return;
      if (!data) {
        setError('Could not start Jupyter server');
        return;
      }

      const host = window.location.hostname || 'localhost';
      const relativePath = resource.startsWith(rootDir) ? resource.slice(rootDir.length).replace(/^\//, '') : resource;
      const url = `http://${host}:${data.port}/doc/tree/${relativePath}`;

      // Poll the notebook URL until it actually serves content
      const deadline = Date.now() + 15000;
      while (Date.now() < deadline && !cancelled) {
        try {
          // Not a call to our server: it probes the jupyter server itself and wants the raw response, so it stays a plain fetch.
          const probe = await fetch(url, { mode: 'no-cors' });
          if (probe.type === 'opaque' || probe.ok) {
            if (!cancelled) setNotebookUrl(url);
            return;
          }
        } catch {
          // not ready yet
        }
        await new Promise(resolve => setTimeout(resolve, 500));
      }
      if (!cancelled) setError('Jupyter notebook took too long to become ready');
    })();

    return () => { cancelled = true; };
  }, [rootDir, resource, themeMode]);

  if (error) return <Typography color="error">{error}</Typography>;
  if (!notebookUrl) return (
    <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%' }}>
      <CircularProgress />
    </Box>
  );

  return (
    <Box sx={{ width: '100%', height: '100%' }}>
      <iframe
        src={notebookUrl}
        style={{ width: '100%', height: '100%', border: 'none' }}
      />
    </Box>
  );
};
