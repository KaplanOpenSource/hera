import { Box, Button, TextField } from '@mui/material';
import { useState } from 'react';
import { API_BASE } from '../shared/constants';
import type { ExecRequest } from '../shared/types';

export const CommandExecutor = () => {
  const [command, setCommand] = useState('');

  const handleExec = async () => {
    const payload: ExecRequest = {
      code: command,
    };
    const r = await fetch(`${API_BASE}/exec`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    const data = await r.json();
    console.log(data);
  };

  return (
    <Box sx={{ mt: 4, display: 'flex', gap: 2 }}>
      <TextField
        label="Command"
        variant="outlined"
        fullWidth
        value={command}
        onChange={(e) => setCommand(e.target.value)}
      />
      <Button variant="contained" onClick={handleExec}>
        Execute
      </Button>
    </Box>
  );
};
