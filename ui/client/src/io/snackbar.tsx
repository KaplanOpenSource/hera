import { forwardRef } from 'react';
import { Box, CircularProgress, IconButton, Paper } from '@mui/material';
import { Close, ErrorOutline, InfoOutlined } from '@mui/icons-material';
import { enqueueSnackbar, closeSnackbar, SnackbarKey } from 'notistack';

type Kind = 'running' | 'error' | 'info';

const ACCENT: Record<Kind, string> = {
  running: '#1976d2',
  error: '#d32f2f',
  info: '#2e7d32',
};

const Toast = forwardRef<HTMLDivElement, {
  kind: Kind,
  message: React.ReactNode,
  onClose: () => void,
}>(({ kind, message, onClose }, ref) => (
  <Paper
    ref={ref}
    elevation={3}
    sx={{
      display: 'flex',
      alignItems: 'center',
      gap: 1,
      pl: 1.25,
      pr: 0.5,
      py: 0.5,
      borderRadius: 999,
      bgcolor: 'rgba(32, 32, 36, 0.92)',
      color: '#fff',
      backdropFilter: 'blur(6px)',
      minWidth: 0,
      maxWidth: 360,
    }}
  >
    {kind === 'running' ? (
      <CircularProgress size={12} thickness={6} sx={{ color: ACCENT.running }} />
    ) : kind === 'info' ? (
      <InfoOutlined sx={{ fontSize: 14, color: ACCENT.info }} />
    ) : (
      <ErrorOutline sx={{ fontSize: 14, color: ACCENT.error }} />
    )}
    <Box sx={{
      fontSize: 12,
      fontWeight: 500,
      lineHeight: 1.3,
      whiteSpace: 'nowrap',
      overflow: 'hidden',
      textOverflow: 'ellipsis',
      flex: 1,
    }}>
      {message}
    </Box>
    <IconButton
      size="small"
      onClick={onClose}
      sx={{ color: 'rgba(255,255,255,0.6)', p: 0.25, '&:hover': { color: '#fff' } }}
    >
      <Close sx={{ fontSize: 14 }} />
    </IconButton>
  </Paper>
));

const toastContent = (kind: Kind) =>
  (key: SnackbarKey, message: React.ReactNode) => (
    <Toast kind={kind} message={message} onClose={() => closeSnackbar(key)} />
  );

export const pushRunning = (label: string): SnackbarKey =>
  enqueueSnackbar(label, {
    persist: true,
    content: toastContent('running'),
  });

export const pushError = (message: string): SnackbarKey =>
  enqueueSnackbar(message, {
    autoHideDuration: 5000,
    content: toastContent('error'),
  });

export const pushInfo = (message: string): SnackbarKey =>
  enqueueSnackbar(message, {
    autoHideDuration: 4000,
    content: toastContent('info'),
  });

export const dismiss = (key: SnackbarKey) => closeSnackbar(key);
