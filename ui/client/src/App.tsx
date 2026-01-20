import { Alert, AppBar, Box, Paper, Stack, Toolbar, Typography } from '@mui/material';
import { useEffect, useState } from 'react';
import { CommandExecutor } from './components/CommandExecutor';
import { DetailsViewPanel } from './components/details/DetailsViewPanel';
import { PageTitle } from './components/header/PageTitle';
import { ProjectChooser } from './components/header/ProjectChooser';
import { ProjectTreeView } from './components/project/ProjectTreeView';
import { FetchProjects } from './io/FetchProjects';
import { useProjectStore } from './stores/useProjectStore';
import { ServerConstantReader } from './stores/useServerConstants';
import { CommitIdShower } from './components/header/CommitIdShower';

export default function App() {
  // const [loading, setLoading] = useState(false);
  // const [error, setError] = useState<string | undefined>(undefined);
  const { getProject } = useProjectStore();
  const [selectedItemsIds, setSelectedItemIds] = useState<string[]>([]);

  const project = getProject();

  useEffect(() => {
    setSelectedItemIds([]);
  }, [project?.name])

  return (<>
    <ServerConstantReader />
    <FetchProjects />
    <Box
      sx={{
        height: '100vh',
        display: 'flex',
        flexDirection: 'column',
      }}
    >
      <AppBar position="static">
        <Toolbar>
          <Stack direction="row" spacing={2}>
            <PageTitle />
            <ProjectChooser />
            <CommitIdShower />
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
          ? (<>
            <Box sx={{ flex: 1, minWidth: 0 }}>
              <Paper sx={{ p: 2, height: '100%', overflow: 'auto' }}>
                <ProjectTreeView
                  project={project}
                  setSelectedItemIds={setSelectedItemIds}
                />
              </Paper>
            </Box>

            <Box sx={{ flex: 1, minWidth: 0 }}>
              <Paper sx={{ p: 2, height: '100%', overflow: 'auto' }}>
                <DetailsViewPanel
                  project={project}
                  showItemId={selectedItemsIds[0]}
                />
              </Paper>
            </Box>
          </>)

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
