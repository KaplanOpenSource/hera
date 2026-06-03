import { Box, Typography } from '@mui/material';

export const PageTitle = () => {
  return (
    <>
      <Box component="img" src="/atom.svg" alt="" sx={{ width: 32, height: 32, display: 'block', alignSelf: 'flex-start' }} />
      <Typography variant="h4" gutterBottom sx={{ marginRight: 2 }}>
        Hera UI
      </Typography>
    </>
  );
};
