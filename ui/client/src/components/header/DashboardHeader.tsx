import { Fullscreen, FullscreenExit, HelpOutline, ViewQuilt } from '@mui/icons-material';
import { AppBar, createTheme, Stack, ThemeProvider, Toolbar } from '@mui/material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { ProjectViewSettingsButton } from '../project/ProjectViewSettingsButton';
import { PageTitle } from './PageTitle';
import { ProjectChooser } from './ProjectChooser';
import { StatusIndicators } from './StatusIndicators';

const headerTheme = createTheme({
  palette: {
    mode: 'dark',
    background: { paper: '#1976d2' },
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
  treeCollapsed,
  setTreeCollapsed,
  onResetLayout,
}: {
  treeCollapsed: boolean,
  setTreeCollapsed: (fn: (c: boolean) => boolean) => void,
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
              title={treeCollapsed ? 'Show sidebar' : 'Hide sidebar'}
              onClick={() => setTreeCollapsed(c => !c)}
              color="inherit"
            >
              {treeCollapsed ? <FullscreenExit /> : <Fullscreen />}
            </ButtonTooltip>
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
            <StatusIndicators />
          </Stack>
        </Toolbar>
      </AppBar>
    </ThemeProvider>
  );
};
