import { Box } from '@mui/material';
import { VersionShower } from './VersionShower';
import { CorsIndicator } from './CorsIndicator';
import { UserIndicator } from './UserIndicator';

export const StatusIndicators = () => {
  return (
    <Box sx={{
      position: 'absolute',
      right: 0,
      top: 0,
      textAlign: 'right',
      lineHeight: 0,
    }}>
      <VersionShower />
      <UserIndicator />
      <CorsIndicator />
    </Box>
  );
};
