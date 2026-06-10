import { Tooltip, Typography } from '@mui/material';
import { fetchPythonClean } from '../../io/fetchPython';
import { useEffect, useState } from 'react';

export const UserIndicator = () => {
  const [username, setUsername] = useState<string | null>(null);

  useEffect(() => {
    (async () => {
      const response = await fetchPythonClean({
        results: ['username'],
        code: 'import getpass; username = getpass.getuser()',
      });
      if (response.data?.username) setUsername(response.data.username);
    })();
  }, []);

  return username && (
    <Tooltip title="User that started the server">
      <Typography
        variant="caption"
        sx={{ fontSize: '10px', cursor: 'default', color: '#7CFC00' }}
      >
        {username}
      </Typography>
    </Tooltip>
  );
};
