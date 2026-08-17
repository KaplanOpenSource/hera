import { HelpOutline, ViewQuilt } from '@mui/icons-material';
import { AppBar, Box, createTheme, Stack, ThemeProvider, Toolbar, Typography } from '@mui/material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { ProjectViewSettingsButton } from '../project/ProjectViewSettingsButton';
import { AddProjectButton } from './AddProjectButton';
import { AutoReloadToggle } from './AutoReloadToggle';
import { CorsIndicator } from './CorsIndicator';
import { PageTitle } from './PageTitle';
import { ProjectChooser } from './ProjectChooser';
import { UserIndicator } from './UserIndicator';
import { VersionShower } from './VersionShower';

const headerTheme = createTheme({
  palette: {
    mode: 'dark',
    primary: { main: '#22d3ee' },
    // Dark navy header bar that blends with the app background (was blue #1976d2).
    background: { paper: '#0b1220' },
  },
  components: {
    MuiInputBase: {
      styleOverrides: {
        input: {
          '&::selection': {
            backgroundColor: 'rgba(255,255,255,0.3)',
            color: '#fff',
          },
        },
      },
    },
  },
});

export const DashboardHeader = ({
  onResetLayout,
}: {
  onResetLayout: () => void,
}) => {
  return (
    <ThemeProvider theme={headerTheme}>
      <AppBar position="static">
        <Toolbar>
          {/* Left: title, project, and the user + CORS indicators. */}
          <Stack direction="row" spacing={1} alignItems="center">
            <PageTitle />
            <ProjectChooser />
            <Stack direction="column" justifyContent="center">
              <UserIndicator />
              <CorsIndicator />
            </Stack>
          </Stack>

          <Box sx={{ flexGrow: 1 }} />

          {/* Right: add project, action buttons, settings, then the version. */}
          <Stack direction="row" spacing={1} alignItems="center">
            <AddProjectButton />
            <ButtonTooltip
              title="Reset panel layout"
              onClick={onResetLayout}
              color="inherit"
            >
              <ViewQuilt />
            </ButtonTooltip>
            <ButtonTooltip
              title="Documentation"
              onClick={() => window.open('https://kaplanopensource.github.io/hera', '_blank')}
              color="inherit"
            >
              <HelpOutline />
            </ButtonTooltip>
            <Stack direction="row" spacing={0.5} alignItems="center">
              <Typography variant="body2">Auto-reload</Typography>
              <AutoReloadToggle />
            </Stack>
            <ProjectViewSettingsButton />
            <VersionShower />
          </Stack>
        </Toolbar>
      </AppBar>
    </ThemeProvider>
  );
};
