import { AppBar, Box, Paper, Stack, Toolbar, Typography } from '@mui/material';
import { useEffect, useState } from 'react';
import { Group as PanelGroup, Panel, Separator as PanelResizeHandle } from 'react-resizable-panels';
import { useParams } from 'react-router-dom';
import { DetailsViewPanel } from './components/details/DetailsViewPanel';
import { PageTitle } from './components/header/PageTitle';
import { ProjectChooser } from './components/header/ProjectChooser';
import { ProjectTreeView } from './components/project/ProjectTreeView';
import { FetchProjects } from './io/FetchProjects';
import { useProjectStore } from './stores/useProjectStore';
import { ServerConstantReader } from './stores/useServerConstants';
import { StatusIndicators } from './components/header/StatusIndicators';

export const Dashboard = () => {
  const { projectName } = useParams<{ projectName: string }>();
  const { getProject } = useProjectStore();
  const [selectedItemsIds, setSelectedItemIds] = useState<string[]>([]);

  const project = getProject();

  useEffect(() => {
    setSelectedItemIds([]);
  }, [project?.name])

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
      <AppBar position="static">
        <Toolbar>
          <Stack direction="row" spacing={2} alignItems="center">
            <PageTitle />
            <ProjectChooser />
            <StatusIndicators />
          </Stack>
        </Toolbar>
      </AppBar>

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
                    setSelectedItemIds={setSelectedItemIds}
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
                <Paper sx={{ p: 2, height: '100%', overflow: 'hidden', display: 'flex', flexDirection: 'column' }}>
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
