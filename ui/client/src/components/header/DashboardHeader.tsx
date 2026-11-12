import { HelpOutline, ViewQuilt } from '@mui/icons-material';
import { AppBar, Box, createTheme, Stack, ThemeProvider, Toolbar, Tooltip } from '@mui/material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { ThemeModeSwitch } from '../../elements/ThemeModeSwitch';
import { ThemeMode, useViewSettingsStore } from '../../stores/useViewSettingsStore';
import { ProjectViewSettingsButton } from '../project/ProjectViewSettingsButton';
import { AddProjectButton } from './AddProjectButton';
import { PageTitle } from './PageTitle';
import { ProjectChooser } from './ProjectChooser';

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
  const { viewSettings, setViewSettings } = useViewSettingsStore();
  return (
    <ThemeProvider theme={headerTheme}>
      <AppBar position="static">
        <Toolbar>
          {/* Left: title, project droplist and its add button. */}
          <Stack direction="row" spacing={1} alignItems="center">
            <PageTitle />
            <ProjectChooser />
            <AddProjectButton />
          </Stack>

          <Box sx={{ flexGrow: 1 }} />

          {/* Right: action buttons, theme switch, then settings. */}
          <Stack direction="row" spacing={1} alignItems="center">
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
            <Tooltip title="Dark mode">
              <span>
                <ThemeModeSwitch
                  mode={viewSettings.themeMode}
                  setMode={(mode: ThemeMode) => setViewSettings({ themeMode: mode })}
                />
              </span>
            </Tooltip>
            <ProjectViewSettingsButton />
          </Stack>
        </Toolbar>
      </AppBar>
    </ThemeProvider>
  );
};
