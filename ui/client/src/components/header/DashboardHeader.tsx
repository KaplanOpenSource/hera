import { HelpOutline, ViewQuilt } from '@mui/icons-material';
import { AppBar, createTheme, Stack, ThemeProvider, Toolbar, Typography } from '@mui/material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { ProjectViewSettingsButton } from '../project/ProjectViewSettingsButton';
import { AutoReloadToggle } from './AutoReloadToggle';
import { PageTitle } from './PageTitle';
import { ProjectChooser } from './ProjectChooser';
import { StatusIndicators } from './StatusIndicators';

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
          <Stack direction="row" spacing={1} alignItems="center">
            <PageTitle />
            <ProjectChooser />
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
            <ProjectViewSettingsButton />
            <Stack direction="row" spacing={0.5} alignItems="center">
              <Typography variant="body2">Auto-reload</Typography>
              <AutoReloadToggle />
            </Stack>
            <StatusIndicators />
          </Stack>
        </Toolbar>
      </AppBar>
    </ThemeProvider>
  );
};
