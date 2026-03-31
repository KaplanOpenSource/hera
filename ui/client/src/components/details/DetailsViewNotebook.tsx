import { Box, CircularProgress, Typography } from '@mui/material';
import { useEffect, useState } from 'react';
import { BASEURL } from '../../shared/baseurl';

export const DetailsViewNotebook = ({
  rootDir,
  notebookName,
}: {
  rootDir: string,
  notebookName: string,
}) => {
  const [jupyterPort, setJupyterPort] = useState<number | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;

    (async () => {
      try {
        const r = await fetch(`${BASEURL}/jupyter/ensure`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ root_dir: rootDir }),
        });
        const data = await r.json();
        if (cancelled) return;

        const host = window.location.hostname || 'localhost';
        const url = `http://${host}:${data.port}/api/status`;
        for (let i = 0; i < 30; i++) {
          try {
            const resp = await fetch(url);
            if (resp.ok) break;
          } catch { /* not ready yet */ }
          if (cancelled) return;
          await new Promise(res => setTimeout(res, 500));
        }

        if (!cancelled) {
          setJupyterPort(data.port);
        }
      } catch {
        if (!cancelled) {
          setError('Could not start Jupyter server');
        }
      }
    })();

    return () => { cancelled = true; };
  }, [rootDir]);

  if (error) return <Typography color="error">{error}</Typography>;
  if (!jupyterPort) return <CircularProgress />;

  const host = window.location.hostname || 'localhost';
  const src = `http://${host}:${jupyterPort}/doc/tree/notebooks/${notebookName}.ipynb`;

  return (
    <Box sx={{ width: '100%', height: '100%' }}>
      <iframe
        src={src}
        style={{ width: '100%', height: '100%', border: 'none' }}
      />
    </Box>
  );
};
