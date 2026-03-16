import { Tooltip, Typography } from '@mui/material';
import { useEffect, useState } from 'react';
import { BASEURL } from '../../shared/baseurl';

export const CorsIndicator = () => {
  const [origins, setOrigins] = useState<string[] | null>(null);

  useEffect(() => {
    (async () => {
      try {
        const url = `${BASEURL}/cors`;
        console.log('CORS fetch url:', url);
        const r = await fetch(url);
        const text = await r.text();
        try {
          const data = JSON.parse(text);
          console.log('CORS response:', data);
          if (data.origins) {
            setOrigins(data.origins);
          }
        } catch {
          console.log('CORS response is not JSON:', text.replace(/\s+/g, ' '));
        }
      } catch (e) {
        console.log('CORS fetch failed:', e);
      }
    })();
  }, []);

  if (!origins) {
    return null;
  }

  return (
    <Tooltip
      title={<>
        Cors origin allowed:
        {origins.map((o) => (
          <div key={o}>{o}</div>
        ))}
      </>}
    >
      <Typography
        variant="caption"
        sx={{
          fontSize: '10px',
          fontWeight: 'bold',
          cursor: 'default',
          color: '#f00',
        }}
      >
        CORS
      </Typography>
    </Tooltip>
  );
};
