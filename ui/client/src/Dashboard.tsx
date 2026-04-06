import { Fullscreen, FullscreenExit } from '@mui/icons-material';
import { AppBar, Box, createTheme, IconButton, Paper, Stack, ThemeProvider, Toolbar, Tooltip, Typography } from '@mui/material';
import { useCallback, useEffect, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import { DetailsViewPanel } from './components/details/DetailsViewPanel';
import { PageTitle } from './components/header/PageTitle';
import { ProjectChooser } from './components/header/ProjectChooser';
import { StatusIndicators } from './components/header/StatusIndicators';
import { ProjectTreeView } from './components/project/ProjectTreeView';
import { SplitWithSidebar } from './elements/SplitWithSidebar';
import { FetchProjects } from './io/FetchProjects';
import { useProjectStore } from './stores/useProjectStore';
import { ServerConstantReader } from './stores/useServerConstants';

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

export const Dashboard = () => {
  const { projectName, docId } = useParams<{ projectName: string; docId: string }>();
  const { getProject } = useProjectStore();
  const navigate = useNavigate();
  const [selectedItemsIds, setSelectedItemIds] = useState<string[]>(
    docId ? [`document_${docId}`] : []
  );
  const [treeCollapsed, setTreeCollapsed] = useState(false);

  const project = getProject();

  useEffect(() => {
    if (docId && project?.documentIds.has(docId)) {
      setSelectedItemIds([`document_${docId}`]);
    } else {
      setSelectedItemIds([]);
    }
  }, [project?.name])

  const handleSetSelectedItemIds = useCallback((ids: string[]) => {
    setSelectedItemIds(ids);
    const selectedId = ids[0];
    const oid = selectedId?.startsWith('document_') ? selectedId.slice('document_'.length) : undefined;
    const basePath = '/' + encodeURIComponent(project?.name ?? '');
    const newPath = oid ? `${basePath}/${oid}` : basePath;
    if (location.pathname !== newPath) {
      navigate(newPath, { replace: true });
    }
  }, [project?.name, navigate]);

  return (<>
    <ServerConstantReader />
    <FetchProjects urlProjectName={projectName} />
    <Box
      sx={{
        height: '100vh',
        display: 'flex',
        flexDirection: 'column',
      }}
    >
      <ThemeProvider theme={headerTheme}>
        <AppBar position="static">
          <Toolbar>
            <Stack direction="row" spacing={1} alignItems="center">
              <PageTitle />
              <ProjectChooser />
              <Tooltip title={treeCollapsed ? 'Show sidebar' : 'Hide sidebar'}>
                <IconButton
                  color="inherit"
                  onClick={() => setTreeCollapsed(c => !c)}
                  size="small"
                >
                  {treeCollapsed ? <FullscreenExit /> : <Fullscreen />}
                </IconButton>
              </Tooltip>
              <StatusIndicators />
            </Stack>
          </Toolbar>
        </AppBar>
      </ThemeProvider>

      <Box
        sx={{
          flex: 1,
          display: 'flex',
          gap: 1,
          minHeight: 0,
        }}
      >
        {project
          ? (
            <SplitWithSidebar
              collapsed={treeCollapsed}
              sidebar={
                <Paper sx={{ p: 2, height: '100%', overflow: 'auto' }}>
                  <ProjectTreeView
                    project={project}
                    selectedItemsIds={selectedItemsIds}
                    setSelectedItemIds={handleSetSelectedItemIds}
                  />
                </Paper>
              }
            >
              <Paper sx={{
                height: '100%',
                overflow: 'hidden',
              }}>
                <DetailsViewPanel
                  project={project}
                  showItemId={selectedItemsIds[0]}
                />
              </Paper>
            </SplitWithSidebar>
          )
          : (
            <Paper sx={{ p: 2, height: '100%', overflow: 'auto', flex: 1, minWidth: 0 }}>
              <Typography>
                No project loaded
              </Typography>
            </Paper>
          )}
      </Box>
    </Box>
  </>)
}
