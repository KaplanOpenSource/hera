import { Box, Tooltip, Typography } from '@mui/material';
import { buildNumber } from '../../buildNumber';
import { fetchPythonClean } from '../../io/fetchPython';
import { useEffect, useState } from 'react';

export const VersionShower = () => {
  const [copied, setCopied] = useState(false);
  const [heraVersion, setHeraVersion] = useState<string | null>(null);

  useEffect(() => {
    (async () => {
      const response = await fetchPythonClean({
        results: ['version'],
        code: 'import hera; version = hera.__version__',
      });
      if (response.data?.version) setHeraVersion(response.data.version);
    })();
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
