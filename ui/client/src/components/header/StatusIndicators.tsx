import { Box } from '@mui/material';
import { CommitIdShower } from './CommitIdShower';
import { CorsIndicator } from './CorsIndicator';

export const StatusIndicators = () => {
  return (
    <Box sx={{ position: 'absolute', right: 0, top: 0, textAlign: 'right', lineHeight: 1 }}>
      <CommitIdShower />
      <CorsIndicator />
    </Box>
  );
};
