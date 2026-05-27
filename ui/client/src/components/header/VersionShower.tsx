import { Box, Tooltip, Typography } from '@mui/material';
import { buildNumber } from '../../buildNumber';
import { BASEURL } from '../../shared/baseurl';
import { useEffect, useState } from 'react';

export const VersionShower = () => {
  const [copied, setCopied] = useState(false);
  const [heraVersion, setHeraVersion] = useState<string | null>(null);

  useEffect(() => {
    const code = 'import hera\nresult = {}\nresult["version"] = hera.__version__';
    fetch(`${BASEURL}/exec`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ code }),
    })
      .then(r => r.json())
      .then(res => {
        if (res.data?.version) setHeraVersion(res.data.version);
      })
      .catch(() => {});
  }, []);

  const versionText = heraVersion
    ? `hera ${heraVersion} | build ${buildNumber}`
    : `build ${buildNumber}`;

  return (
    <Box>
      <Tooltip
        title={<>
          {heraVersion && <>Hera v{heraVersion}<br /></>}
          UI build {buildNumber}
          <br />
          {copied ? 'Copied!' : 'click to copy'}
        </>}
      >
        <Typography
          variant="caption"
          onClick={() => {
            setCopied(true);
            navigator.clipboard.writeText(versionText);
            setTimeout(() => {
              setCopied(false);
            }, 2000);
          }}
          sx={{ fontSize: '10px' }}
        >
          {versionText}
        </Typography>
      </Tooltip>
    </Box>
  );
};
