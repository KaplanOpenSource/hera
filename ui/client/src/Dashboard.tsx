import { AppBar, Box, Paper, Stack, Toolbar, Typography, createTheme, ThemeProvider } from '@mui/material';
import { useCallback, useEffect, useState } from 'react';
import { Group as PanelGroup, Panel, Separator as PanelResizeHandle } from 'react-resizable-panels';
import { useNavigate, useParams } from 'react-router-dom';
import { DetailsViewPanel } from './components/details/DetailsViewPanel';
import { idFromNotebookId } from './shared/idDocId';
import { PageTitle } from './components/header/PageTitle';
import { ProjectChooser } from './components/header/ProjectChooser';
import { ProjectTreeView } from './components/project/ProjectTreeView';
import { FetchProjects } from './io/FetchProjects';
import { useProjectStore } from './stores/useProjectStore';
import { ServerConstantReader } from './stores/useServerConstants';
import { StatusIndicators } from './components/header/StatusIndicators';

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
            <Stack direction="row" spacing={2} alignItems="center">
              <PageTitle />
              <ProjectChooser />
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
          minHeight: 0, // important for scroll behavior
        }}
      >
        {project
          ? (
            <PanelGroup orientation="horizontal">
              <Panel defaultSize={50} minSize={20}>
                <Paper sx={{ p: 2, height: '100%', overflow: 'auto' }}>
                  <ProjectTreeView
                    project={project}
                    selectedItemsIds={selectedItemsIds}
                    setSelectedItemIds={handleSetSelectedItemIds}
                  />
                </Paper>
              </Panel>

              <PanelResizeHandle
                style={{
                  width: 4,
                  cursor: 'col-resize',
                  backgroundColor: '#e0e0e0',
                  outline: 'none',
                }}
              />

              <Panel defaultSize={50} minSize={20}>
                <Paper sx={{
                  p: idFromNotebookId(selectedItemsIds[0]) ? 0 : 2,
                  height: '100%',
                  overflow: 'hidden',
                }}>
                  <DetailsViewPanel
                    project={project}
                    showItemId={selectedItemsIds[0]}
                  />
                </Paper>
              </Panel>
            </PanelGroup>
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
