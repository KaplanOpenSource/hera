import { Box, Tooltip, Typography } from '@mui/material';
import { commitId } from '../../commitId';
import { useState } from 'react';

export const CommitIdShower = () => {
  const [copied, setCopied] = useState(false);
  return (
    <Box sx={{ position: 'absolute', right: 0, top: 0 }}>
      <Tooltip
        title={<>
          UI build commit id
          <br />
          {commitId}
          <br />
          {copied ? 'Copied!' : 'click to copy'}
        </>}
      >
        <Typography
          variant="caption"
          onClick={() => {
            setCopied(true);
            navigator.clipboard.writeText(commitId);
            setTimeout(() => {
              setCopied(false);
            }, 2000);
          }}
          sx={{fontSize: '10px'}}
        >
          {commitId.substring(0, 6)}
        </Typography>
      </Tooltip>
    </Box>
  );
};
