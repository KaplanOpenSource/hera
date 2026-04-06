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
    setJupyterPort(null);
    setError(null);

    (async () => {
      try {
        const r = await fetch(`${BASEURL}/jupyter/ensure`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ root_dir: rootDir }),
        });
        const data = await r.json();
        if (cancelled) return;

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
  }, [rootDir, notebookName]);

  if (error) return <Typography color="error">{error}</Typography>;
  if (!jupyterPort) return (
    <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%' }}>
      <CircularProgress />
    </Box>
  );

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
