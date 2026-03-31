import { Box, CircularProgress, Typography } from '@mui/material';
import { useEffect, useState } from 'react';
import { BASEURL } from '../../shared/baseurl';

export const DetailsViewNotebook = ({
  name,
}: {
  name: string,
}) => {
  const [jupyterPort, setJupyterPort] = useState<number | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    (async () => {
      try {
        const r = await fetch(`${BASEURL}/jupyter`);
        const data = await r.json();
        setJupyterPort(data.port);
      } catch {
        setError('Could not reach server');
      }
    })();
  }, []);

  if (error) return <Typography color="error">{error}</Typography>;
  if (!jupyterPort) return <CircularProgress />;

  const host = window.location.hostname || 'localhost';
  const src = `http://${host}:${jupyterPort}/lab`;

  return (
    <Box sx={{ width: '100%', height: '100%' }}>
      <iframe
        src={src}
        style={{ width: '100%', height: '100%', border: 'none' }}
      />
    </Box>
  );
};
