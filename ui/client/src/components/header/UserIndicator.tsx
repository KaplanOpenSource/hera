import { Tooltip, Typography } from '@mui/material';
import { useEffect } from 'react';
import { useServerUserStore } from '../../stores/useServerUserStore';

export const UserIndicator = () => {
  const { username, inDocker, loadServerUser } = useServerUserStore();

  useEffect(() => {
    loadServerUser();
  }, [loadServerUser]);

  return username && (
    <Tooltip title="User that started the server">
      <Typography
        variant="caption"
        sx={{ fontSize: '10px', cursor: 'default', color: '#4ade80' }}
      >
        {username}
        {inDocker && (
          <span style={{ color: 'rgba(74, 222, 128, 0.5)' }}> (docker)</span>
        )}
      </Typography>
    </Tooltip>
  );
};
