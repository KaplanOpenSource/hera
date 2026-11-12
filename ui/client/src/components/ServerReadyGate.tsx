import { Box, CircularProgress, Typography } from '@mui/material';
import { ReactNode, useEffect, useState } from 'react';
import { useLocation, useNavigate } from 'react-router-dom';
import { fetchInit } from '../io/fetchInit';
import { projectPath } from '../shared/projectPath';

const READY_POLL_MS = 1000;

// The project name a path names, if it names one.
const projectFromPath = (pathname: string) => {
  const [, first] = pathname.split('/');
  return first ? decodeURIComponent(first) : undefined;
};

// Spinner until the server warmed hera; polls /init, which answers with
// everything the UI needs to start, then mounts the app with the data in place.
export const ServerReadyGate = ({
  children,
}: {
  children: ReactNode,
}) => {
  const [ready, setReady] = useState(false);
  const navigate = useNavigate();
  // The router's path, not window.location: under a router that is not the
  // browser history the window shows none of this, and asking /init for the
  // wrong project made it answer with another project and navigate to it.
  const { pathname } = useLocation();

  useEffect(() => {
    if (ready) {
      return;
    }
    let cancelled = false;
    let timer: ReturnType<typeof setTimeout>;

    console.log('[ServerReadyGate] waiting for server to warm up…');
    const poll = async () => {
      const wanted = projectFromPath(pathname);
      let data = null;
      try {
        data = await fetchInit(wanted);
      } catch (e) {
        // Not answering yet — log and keep polling.
        console.error('[ServerReadyGate] init failed:', e);
      }
      if (cancelled) {
        return;
      }
      if (data?.ready) {
        if (data.problems) {
          console.error('[ServerReadyGate] init problems:', data.problems);
        }
        console.log('[ServerReadyGate] server ready — mounting app');
        if (data.project && data.project.name !== wanted) {
          navigate(projectPath(data.project.name), { replace: true });
        }
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
  }, [ready, navigate, pathname]);

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
