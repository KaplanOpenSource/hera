import { Box, Tooltip, Typography } from '@mui/material';
import { buildNumber } from '../../buildNumber';
import { useState } from 'react';

export const CommitIdShower = () => {
  const [copied, setCopied] = useState(false);
  return (
    <Box>
      <Tooltip
        title={<>
          UI build number
          <br />
          {buildNumber}
          <br />
          {copied ? 'Copied!' : 'click to copy'}
        </>}
      >
        <Typography
          variant="caption"
          onClick={() => {
            setCopied(true);
            navigator.clipboard.writeText(buildNumber);
            setTimeout(() => {
              setCopied(false);
            }, 2000);
          }}
          sx={{ fontSize: '10px' }}
        >
          {buildNumber}
        </Typography>
      </Tooltip>
    </Box>
  );
};
