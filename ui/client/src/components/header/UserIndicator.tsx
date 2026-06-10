import { Tooltip, Typography } from '@mui/material';
import { fetchPythonClean } from '../../io/fetchPython';
import { useEffect, useState } from 'react';

export const UserIndicator = () => {
  const [username, setUsername] = useState<string | null>(null);
  const [inDocker, setInDocker] = useState(false);

  useEffect(() => {
    (async () => {
      const response = await fetchPythonClean({
        results: ['username', 'inDocker'],
        code: "import getpass, os; username = getpass.getuser(); inDocker = os.path.exists('/.dockerenv')",
      });
      if (response.data?.username) setUsername(response.data.username);
      setInDocker(Boolean(response.data?.inDocker));
    })();
  }, []);

  return username && (
    <Tooltip title="User that started the server">
      <Typography
        variant="caption"
        sx={{ fontSize: '10px', cursor: 'default', color: '#7CFC00' }}
      >
        {username}
        {inDocker && (
          <span style={{ color: 'rgba(124, 252, 0, 0.5)' }}> (docker)</span>
        )}
      </Typography>
    </Tooltip>
  );
};
