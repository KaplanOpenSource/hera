import { Box, Button, Stack, TextField, Tooltip } from '@mui/material';
import { useState } from 'react';
import { execPython } from '../io/execPython';

export const CommandExecutor = () => {
  const [command, setCommand] = useState('');

  const handleExec = async () => {
    // const { data } = 
    await execPython(command)
    // console.log(data);
  };

  return (
    <Tooltip
      title={<span>Executes a python command<br />The variable named result will be printed to console</span>}
    >
      <Stack direction={'row'} spacing={1} sx={{ bottom: 10, right: 10, left: 10, position: 'absolute' }}>
        <TextField
          // size='small'
          label="Command"
          variant="outlined"
          fullWidth
          value={command}
          onChange={(e) => setCommand(e.target.value)}
        />
        <Button
          size='small'
          variant="contained"
          sx={{ textTransform: 'none' }}
          onClick={handleExec}
        >
          Execute Python
        </Button>
      </Stack>
    </Tooltip>
  );
};
