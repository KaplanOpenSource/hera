import { Box } from '@mui/material';

const JUPYTER_PORT = 8888;

export const DetailsViewNotebook = ({
  name,
}: {
  name: string,
}) => {
  const host = window.location.hostname || 'localhost';
  const src = `http://${host}:${JUPYTER_PORT}/notebooks/${encodeURIComponent(name)}.ipynb`;

  return (
    <Box sx={{ width: '100%', height: '100%' }}>
      <iframe
        src={src}
        style={{ width: '100%', height: '100%', border: 'none' }}
      />
    </Box>
  );
};
