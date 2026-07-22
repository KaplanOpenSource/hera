import { Box, CircularProgress, Typography } from '@mui/material';
import { useEffect, useState } from 'react';
import { BASEURL } from '../../shared/baseurl';
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
      try {
        const r = await fetch(`${BASEURL}/jupyter/ensure`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ root_dir: rootDir, dark: themeMode === ThemeMode.Dark }),
        });
        const data = await r.json();
        if (cancelled) return;

        const host = window.location.hostname || 'localhost';
        const relativePath = resource.startsWith(rootDir) ? resource.slice(rootDir.length).replace(/^\//, '') : resource;
        const url = `http://${host}:${data.port}/doc/tree/${relativePath}`;

        // Poll the notebook URL until it actually serves content
        const deadline = Date.now() + 15000;
        while (Date.now() < deadline && !cancelled) {
          try {
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
      } catch {
        if (!cancelled) setError('Could not start Jupyter server');
      }
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
