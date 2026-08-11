import { Box, CircularProgress, Typography } from '@mui/material';
import { ReactNode, useEffect, useState } from 'react';
import { BASEURL } from '../shared/baseurl';

const READY_POLL_MS = 1000;

// Spinner until the server warmed hera; polls /ready, then mounts the app.
export const ServerReadyGate = ({
  children,
}: {
  children: ReactNode,
}) => {
  const [ready, setReady] = useState(false);

  useEffect(() => {
    if (ready) {
      return;
    }
    let cancelled = false;
    let timer: ReturnType<typeof setTimeout>;

    console.log('[ServerReadyGate] waiting for server to warm up…');
    const poll = async () => {
      let isReady = false;
      try {
        const r = await fetch(`${BASEURL}/ready`);
        const contentType = r.headers.get('content-type');
        // Old server (no /ready route) answers with the SPA's HTML — treat as not ready, keep polling.
        if (r.ok && contentType?.includes('application/json')) {
          const data = await r.json();
          isReady = Boolean(data.ready);
        }
        console.log(`[ServerReadyGate] poll: ${r.status}, ready=${isReady}`);
      } catch (e) {
        // Not answering yet — log and keep polling.
        console.error('[ServerReadyGate] poll failed:', e);
      }
      if (cancelled) {
        return;
      }
      if (isReady) {
        console.log('[ServerReadyGate] server ready — mounting app');
        setReady(true);
      } else {
        timer = setTimeout(poll, READY_POLL_MS);
      }
    };
    poll();

    return () => {
      cancelled = true;
      clearTimeout(timer);
    };
  }, [ready]);

  return ready
    ? <>{children}</>
    : (
      <Box
        sx={{
          height: '100vh',
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          justifyContent: 'center',
          gap: 2,
        }}
      >
        <CircularProgress />
        <Typography>Starting Hera…</Typography>
      </Box>
    );
};
