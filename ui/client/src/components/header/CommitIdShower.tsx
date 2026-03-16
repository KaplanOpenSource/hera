import { Box, Tooltip, Typography } from '@mui/material';
import { commitId } from '../../commitId';
import { useState } from 'react';

export const CommitIdShower = () => {
  const [copied, setCopied] = useState(false);
  const dev = import.meta.env.MODE === 'development';
  const id = dev ? 'dev' : commitId;
  return (
    <Box>
      <Tooltip
        title={dev
          ? <>
            UI is running in development
          </>
          : <>
            UI build commit id
            <br />
            {id}
            <br />
            {copied ? 'Copied!' : 'click to copy'}
          </>}
      >
        <Typography
          variant="caption"
          onClick={() => {
            if (!dev) {
              setCopied(true);
              navigator.clipboard.writeText(id);
              setTimeout(() => {
                setCopied(false);
              }, 2000);
            }
          }}
          sx={{ fontSize: '10px' }}
        >
          {id.substring(0, 6)}
        </Typography>
      </Tooltip>
    </Box>
  );
};
