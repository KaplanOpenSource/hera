import { Tooltip, Typography } from '@mui/material';
import { useEffect, useState } from 'react';
import { serverApi } from '../../api/api';

export const CorsIndicator = () => {
  const [origins, setOrigins] = useState<string[] | null>(null);

  useEffect(() => {
    (async () => {
      const data = await serverApi.cors();
      if (data?.origins) {
        setOrigins(data.origins);
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
        {origins.includes('*') ? 'CORS ALL' : 'CORS'}
      </Typography>
    </Tooltip>
  );
};
